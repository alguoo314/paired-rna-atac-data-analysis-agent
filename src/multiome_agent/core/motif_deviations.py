"""Fixed-core chromVAR-style motif deviation scores.

Uses `pychromvar` for the actual chromVAR algorithm (motif matching via
MOODS, background-peak sampling, bias-corrected per-cell deviation z-scores)
-- but NOT its own `get_genome()`/`add_peak_seq()` genome download, which
would pull a second, differently-formatted copy of hg38 we already have as
`data/raw/reference/hg38.2bit`. Instead we populate `adata.uns["peak_seq"]`
ourselves via `py2bit` reads against that file, then hand off to
pychromvar's `add_gc_bias`/`match_motif`/`get_bg_peaks`/`compute_deviations`,
which only need `peak_seq` to already be populated -- verified by reading
their source.

Motifs: the full JASPAR 2024 CORE vertebrates non-redundant set (879
motifs, ~280KB) -- this is already the curated, redundancy-collapsed set
JASPAR itself recommends for scans like this, so no further hand-picking on
top of it. Fetched once and cached under `data/raw/reference/` (gitignored).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import py2bit
import pychromvar as pc
import urllib.request
from anndata import AnnData
from Bio import motifs as biomotifs
from mudata import MuData

from multiome_agent.config import REPO_ROOT
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

GENOME_2BIT = REPO_ROOT / "data" / "raw" / "reference" / "hg38.2bit"
JASPAR_URL = (
    "https://jaspar.elixir.no/download/data/2024/CORE/"
    "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
)
JASPAR_CACHE = REPO_ROOT / "data" / "raw" / "reference" / "jaspar2024_core_vertebrates_nr.jaspar"


def _fetch_jaspar_motifs() -> list:
    if not JASPAR_CACHE.exists():
        logger.info("Downloading JASPAR 2024 CORE vertebrates non-redundant motifs to %s", JASPAR_CACHE)
        JASPAR_CACHE.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(JASPAR_URL, JASPAR_CACHE)
    with open(JASPAR_CACHE) as f:
        motif_list = list(biomotifs.parse(f, "jaspar"))
    logger.info("Loaded %d JASPAR motifs", len(motif_list))
    return motif_list


def _add_peak_seq_from_2bit(atac: AnnData) -> AnnData:
    """Populate `atac.uns["peak_seq"]` from `hg38.2bit`, dropping peaks whose
    chrom isn't in the 2bit file (unplaced/alt scaffolds use a different
    naming convention there than in the Cell Ranger ARC reference; a small,
    documented loss of coverage rather than a naming-reconciliation project).
    """
    tb = py2bit.open(str(GENOME_2BIT))
    known_chroms = set(tb.chroms())
    keep = atac.var["chrom"].isin(known_chroms)
    dropped = int((~keep).sum())
    if dropped:
        logger.info("Dropping %d/%d peaks not present in hg38.2bit's chrom set (alt/unplaced scaffolds)", dropped, atac.n_obs)
    atac = atac[:, keep].copy()

    seqs = [
        tb.sequence(row.chrom, int(row.start), int(row.end)).upper()
        for row in atac.var.itertuples()
    ]
    tb.close()
    atac.uns["peak_seq"] = seqs
    return atac


def compute_motif_deviations(mdata: MuData) -> AnnData:
    """Compute per-cell chromVAR-style motif deviation z-scores.

    Returns an AnnData (cells x motifs) of deviation scores, also stashed at
    `mdata.mod["atac"].uns["chromvar_motif_names"]` (motif name list) and
    `mdata.mod["atac"].obsm["chromvar_deviations"]` (same values as a
    DataFrame, motif names as columns, for easy citation by name).
    """
    atac = mdata.mod["atac"].copy()
    atac.X = atac.layers["counts"] if "counts" in atac.layers else atac.X
    # `compute_deviations` divides `adata.X` in place (`a /= np.sum(count)`),
    # which numpy rejects when `X` is an unsigned-integer dtype (float result
    # can't safely cast back into uint) -- raised only on data whose raw
    # counts are uint (the shareseq-multi-cell-lines dataset's ATAC `.X`; the tenx-cell-ranger
    # loader's `.X` happens to already be float32). Cast defensively so this
    # works regardless of the source dtype, not just the one dataset that
    # happened to already be float.
    atac.X = atac.X.astype(np.float64)

    atac = _add_peak_seq_from_2bit(atac)
    pc.add_gc_bias(atac)

    motif_list = _fetch_jaspar_motifs()
    pc.match_motif(atac, motif_list, background="even")

    logger.info("Sampling GC/accessibility-matched background peaks...")
    pc.get_bg_peaks(atac)

    logger.info("Computing chromVAR deviations for %d motifs x %d cells...", len(motif_list), atac.n_obs)
    dev = pc.compute_deviations(atac)

    dev_df = pd.DataFrame(np.asarray(dev.X), index=dev.obs_names, columns=dev.var_names)
    dev_df = dev_df.reindex(mdata.mod["atac"].obs_names)

    mdata.mod["atac"].uns["chromvar_motif_names"] = list(dev.var_names)
    mdata.mod["atac"].obsm["chromvar_deviations"] = dev_df

    logger.info("Motif deviations complete: %d motifs retained", dev.n_vars)
    return dev


def _motif_match_rank(query: str, column: str) -> int:
    """Lower is a better match. JASPAR chromvar_deviations columns look like
    "MA0080.7.Spi1" or, for heterodimers, "MA0019.2.Ddit3::Cebpa" -- the real
    TF-name segment is everything after the matrix ID, further split on
    "::" for composites. 0 = `query` exactly equals that segment alone; 1 =
    `query` is exactly one factor of a composite; 2 = only a raw substring
    match (the previous, sole matching criterion) -- kept as a last-resort
    fallback, but never preferred over an exact token match.

    Added after a real bug an adversarial "judge" agent caught in this
    session: searching for "EBF1" substring-matched "MA0595.1.SREBF1" (which
    contains "ebf1") ahead of the real, exact "MA0154.5.EBF1" match, because
    the old code took `matches[0]` with no ranking at all -- whichever
    happened to sort first. Verified against real chromvar_deviations
    columns, not assumed: this dataset has both.
    """
    name_part = column.split(".", 2)[-1] if column.count(".") >= 2 else column
    tokens = [t.lower() for t in re.split(r"::|_", name_part)]
    query_l = query.lower()
    if tokens == [query_l]:
        return 0
    if query_l in tokens:
        return 1
    return 2


def tf_expression_motif_correlation(mdata: MuData, gene: str, motif_name_contains: str) -> dict:
    """Spearman correlation between a TF's own RNA expression and its own
    motif's chromVAR deviation, across cells -- a lightweight sanity check
    for CLAUDE.md's "TF expression vs. TF motif accessibility" principle.

    Spearman (rank-based) is invariant to log1p, so it doesn't matter whether
    `rna.X` is raw or log-normalized at this point -- but it DOES matter
    whether per-cell library-size normalization has already been applied
    (that reweights cells differently, so it changes rank order); call this
    after `cluster_rna`'s `normalize_total`, not before.

    A moderate positive correlation is the expected result, not a strict
    pass/fail bar -- CLAUDE.md is explicit that RNA-ATAC relationships are
    noisy per gene. This is a confidence check, not a gate.
    """
    from scipy.stats import spearmanr

    rna = mdata.mod["rna"]
    if gene not in rna.var_names:
        raise ValueError(f"{gene!r} not found in RNA var_names")
    expr = np.asarray(rna[:, gene].X.todense()).ravel() if hasattr(rna[:, gene].X, "todense") else np.asarray(rna[:, gene].X).ravel()

    dev_df = mdata.mod["atac"].obsm["chromvar_deviations"]
    matches = [c for c in dev_df.columns if motif_name_contains.lower() in c.lower()]
    if not matches:
        raise ValueError(f"No motif matching {motif_name_contains!r} among chromvar_deviations columns")
    motif_name = min(matches, key=lambda c: _motif_match_rank(motif_name_contains, c))
    deviation = dev_df[motif_name].reindex(rna.obs_names).to_numpy()

    rho, pval = spearmanr(expr, deviation)
    logger.info(
        "TF-motif sanity check: %s expression vs %s deviation: rho=%.3f (p=%.2e)",
        gene, motif_name, rho, pval,
    )
    result = {"gene": gene, "motif": motif_name, "spearman_rho": float(rho), "pvalue": float(pval)}
    mdata.uns["tf_motif_sanity_check"] = result
    return result
