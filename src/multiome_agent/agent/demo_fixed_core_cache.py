"""A small, deliberately-trimmed fixed-core object committed directly to
the repo (`demo_cache/agent_demo_core.h5mu`, ~46MB) so `make demo` is
genuinely fast on a fresh clone -- unlike the full report-generation
pipeline (`agent/fixed_core_cache.py`'s `get_agent_fixed_core`: 11,909
cells, ~2.6GB, never committed, since GitHub hard-rejects any git-tracked
file over 100MB and this project's own real full-scale rebuild takes
60-90+ minutes).

Committable at all only because it's trimmed to exactly what the 13 agent
tools (`agent/loop.py`'s `TOOLS`) read at query time -- confirmed by
grepping every tool's implementation, not guessed: RNA log-normalized
expression + QC/cluster `.obs` columns + `rank_genes_groups`; ATAC
QC/cluster `.obs` columns + gene-activity + chromVAR-deviations matrices +
(added for `peak_to_gene_links`/`regulon_inference`, Phase 2) real peak
coordinates (`var["chrom"/"start"/"end"]`), the real per-cell peak-count
matrix (`layers["counts"]`), and the peak x motif match matrix
(`varm["motif_match"]`/`uns["motif_match_names"]`). Dropped: PCA/LSI
loadings and descriptive peak/gene annotation columns beyond
chrom/start/end (clustering scratch needed only to BUILD the cluster
labels/DE tables that are kept, never read again afterward). This shrunk a
naive n=150 rebuild from 96.8MB to 46.25MB before the Phase 2 additions
above added the real ATAC counts matrix back (needed for genuine
peak-to-gene/regulon correlations, not just gene-activity/chromVAR, which
only ever consumed it as an intermediate). Also strips `atac.uns["files"]`,
which otherwise bakes this machine's absolute local fragments-file path
into the committed binary -- a real leak caught before committing, not a
hypothetical one.

A real gap caught after the fact, not before: this trim function predated
`peak_to_gene_links`/`regulon_inference` entirely, and originally stripped
every one of the four fields they need (chrom/start/end wiped to zero
columns, varm cleared, motif_match_names excluded from the uns keep-list,
counts layer dropped and `.X` zeroed) -- so the committed demo cache could
offer those two tools in its tool list while every real call to either
failed. Fixed by adding the four fields above; the already-committed demo
cache was separately repaired via `scripts/backfill_demo_regulon_fields.py`
(reindexing from the already-built, already-motif_match-persisted real
agent fixed-core cache's matching 150 cell barcodes / matching peak names,
rather than paying this module's own ~17-20 min from-scratch rebuild cost
again) -- see that script's docstring for the one-time repair, and
PROGRESS_phase2.md for the full story.

Illustrative only, at a small subsample deliberately not meant to support
rigorous claims -- see `reports/examples/` for the real, full-scale
analyses this project's flagship reports are built from.
"""

from __future__ import annotations

import scipy.sparse as sp
from mudata import MuData, read_h5mu, write_h5mu

from multiome_agent.config import REPO_ROOT
from multiome_agent.core.motif_deviations import motif_name_tokens
from multiome_agent.core.pipeline import run_fixed_core
from multiome_agent.data.loader import RAW_DIR, load_pbmc_multiome
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

DEMO_N_CELLS = 150
DEMO_SEED = 0
FRAGMENTS_FILE = RAW_DIR / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
DEMO_CACHE_PATH = REPO_ROOT / "demo_cache" / "agent_demo_core.h5mu"


def _trim_to_tool_readable_fields(mdata: MuData) -> MuData:
    """Keeps only what `agent/loop.py`'s `TOOLS` actually read at query
    time (see module docstring) -- everything else is clustering/embedding
    scratch or purely descriptive annotation that's dead weight once
    cluster labels and DE/gene-activity/chromVAR results are already baked
    in, plus one real path leak (`atac.uns["files"]`).

    `atac.var` is kept down to just `chrom`/`start`/`end` (not the full
    descriptive annotation) -- the minimum `peak_to_gene_links`/
    `regulon_inference` need to locate peaks genomically, dropping
    `gene_ids`/`feature_types`/`genome`/`interval`/`nearest_genes`/
    `distance`/`peak_types`/`primary_gene`/`n_cells`, which no tool reads.
    `atac.layers["counts"]` (real per-cell peak accessibility) is kept
    whole -- both tools need real counts to correlate against. `atac.varm
    ["motif_match"]`/`atac.uns["motif_match_names"]` are kept only for
    motifs whose real TF-name token (`core.motif_deviations.motif_name_tokens`)
    matches a gene actually present in this RNA panel -- `regulon_inference`
    rejects any OTHER TF at its very first gate (`tf_gene not in
    rna.var_names`) regardless of whether its motif column exists, so
    dropping those columns changes zero real behavior while cutting the
    matrix from 879 to ~400-something columns (JASPAR's non-human/
    non-expressed/composite-only entries are the majority of what's
    dropped) -- the dominant single contributor to file size at this
    matrix's real ~20%+ density (where sparse CSR storage costs MORE than
    dense uint8, verified by measuring the actual h5 layout, not assumed).
    `.X` stays zeroed since both tools prefer `layers["counts"]` and never
    fall through to it once that layer exists.
    """
    rna, atac = mdata.mod["rna"], mdata.mod["atac"]

    rna.var = rna.var[[]]
    rna.varm.clear()
    rna.obsm.clear()
    for key in list(rna.uns.keys()):
        if key != "rank_genes_groups":
            del rna.uns[key]
    rna.layers.pop("counts", None)

    atac.var = atac.var[["chrom", "start", "end"]]
    for key in list(atac.varm.keys()):
        if key != "motif_match":
            del atac.varm[key]
    for key in ("X_lsi", "X_lsi_no_depth"):
        atac.obsm.pop(key, None)
    for key in list(atac.uns.keys()):
        if key not in ("gene_activity_genes", "chromvar_motif_names", "motif_match_names"):
            del atac.uns[key]  # drops uns["files"]'s local fragments path too
    atac.X = sp.csr_matrix((atac.n_obs, atac.n_vars))  # unused by any tool; kept shape-valid, zero storage cost

    if "motif_match" in atac.varm and "motif_match_names" in atac.uns:
        rna_genes = {g.upper() for g in rna.var_names}  # index, not a column -- survives var[[]] above
        motif_names = list(atac.uns["motif_match_names"])
        keep_idx = [
            i for i, name in enumerate(motif_names)
            if any(t.upper() in rna_genes for t in motif_name_tokens(name))
        ]
        atac.varm["motif_match"] = atac.varm["motif_match"][:, keep_idx]
        atac.uns["motif_match_names"] = [motif_names[i] for i in keep_idx]

    mdata.update()
    return mdata


def build_demo_cache(n_cells: int = DEMO_N_CELLS, seed: int = DEMO_SEED) -> MuData:
    """Builds the small, trimmed, committable demo cache from scratch and
    writes it to `DEMO_CACHE_PATH`. Real work -- still pays the ~17-20 min
    fixed fragments-import cost (see `core/pipeline.py`) -- run once by a
    maintainer and committed, never run by end users. See
    `scripts/build_demo_cache.py` for the CLI entry point.
    """
    mdata = load_pbmc_multiome(n_cells=n_cells, seed=seed, force_reload=False)
    mdata = run_fixed_core(mdata, FRAGMENTS_FILE)
    mdata = _trim_to_tool_readable_fields(mdata)
    DEMO_CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    write_h5mu(DEMO_CACHE_PATH, mdata)
    logger.info("Wrote demo cache to %s", DEMO_CACHE_PATH)
    return mdata


def get_demo_fixed_core() -> MuData:
    """Loads the committed demo cache. Errors clearly if missing rather
    than silently rebuilding -- rebuilding would defeat the entire point
    (`make demo` must never trigger the ~17-20 min pipeline on an end
    user's machine)."""
    if not DEMO_CACHE_PATH.exists():
        raise FileNotFoundError(
            f"{DEMO_CACHE_PATH} not found. This is meant to be committed to the repo -- "
            "if you're intentionally regenerating it, run scripts/build_demo_cache.py."
        )
    return read_h5mu(DEMO_CACHE_PATH)
