"""Peak-to-gene links (cis-co-accessibility, Cicero/ArchR-style): does a
DISTAL peak's accessibility track a gene's own RNA expression, across
cells? This is a real test for CLAUDE.md's "distal-enhancer regulation"
explanation (core principle 2) -- previously something the agent was only
*allowed* to invoke by name when RNA/ATAC disagreed, never something a tool
actually checked.

Deliberately scoped to ONE gene's own expression on one side of the
correlation, and ONE candidate peak's own accessibility on the other --
never a transcription factor's expression standing in for the gene (that's
a different question this module does not answer). "Distal" is defined
relative to `core.gene_activity`'s own promoter-extended gene-activity
window: a peak already inside that window is what gene activity already
measures, so it's explicitly excluded here -- a real distal-peak
accessibility score should tell you something gene activity doesn't.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from mudata import MuData
from scipy.stats import false_discovery_control, spearmanr

from multiome_agent.core.gene_activity import load_protein_coding_gene_coords, promoter_extend_gene_windows
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

DEFAULT_WINDOW_BP = 500_000
DEFAULT_PROXIMAL_UPSTREAM_BP = 2000  # matches gene_activity's own promoter-extension default
# Minimum |rho| to call a link "significant" alongside the BH q-value gate --
# q-value alone isn't a safe bar at real single-cell sample sizes (thousands
# of cells make even a trivial rho reach q<<0.05; see _correlate_candidate_peaks'
# docstring for the real PBMC numbers that motivated this).
DEFAULT_MIN_ABS_RHO = 0.2


def _select_distal_candidate_peaks(
    peaks: pd.DataFrame, gene_chrom: str, gene_start: int, gene_end: int, gene_strand: str,
    window_bp: int = DEFAULT_WINDOW_BP, proximal_upstream_bp: int = DEFAULT_PROXIMAL_UPSTREAM_BP,
) -> pd.DataFrame:
    """Peaks within `window_bp` of the gene body, on the gene's own
    chromosome, EXCLUDING anything inside the same promoter-extended window
    `core.gene_activity._load_protein_coding_gene_windows` would already
    count toward gene activity -- the point is to test accessibility gene
    activity does NOT already capture.

    `peaks` must have `chrom`/`start`/`end` columns in the SAME ROW ORDER
    as the ATAC matrix's columns (e.g. `atac.var[["chrom", "start", "end"]]`
    passed through unmodified) -- row position, not the pandas index label,
    is what lets the caller slice the right accessibility column per
    surviving candidate. Returns a new frame (own `_peak_pos` column
    recording that original row position, plus `distance_bp`) restricted to
    chrom/window/proximal-exclusion matches -- empty if none.

    Kept as a pure function of plain arrays/frames (no mdata, no genome
    file) so the overlap/exclusion logic is unit-testable on small
    synthetic coordinate tables.
    """
    on_chrom = peaks.reset_index(drop=True)
    on_chrom = on_chrom[on_chrom["chrom"] == gene_chrom]
    if on_chrom.empty:
        return on_chrom.assign(_peak_pos=pd.Series(dtype="int64"), distance_bp=pd.Series(dtype="int64"))

    prox_start = gene_start - proximal_upstream_bp if gene_strand == "+" else gene_start
    prox_end = gene_end if gene_strand == "+" else gene_end + proximal_upstream_bp
    prox_start = max(prox_start, 0)

    search_start = max(gene_start - window_bp, 0)
    search_end = gene_end + window_bp

    p_start, p_end = on_chrom["start"].to_numpy(), on_chrom["end"].to_numpy()
    in_search_window = (p_start < search_end) & (p_end > search_start)
    in_proximal_zone = (p_start < prox_end) & (p_end > prox_start)
    keep = in_search_window & ~in_proximal_zone

    distance_bp = np.where(
        p_start > gene_end, p_start - gene_end,
        np.where(p_end < gene_start, gene_start - p_end, 0),
    )
    result = on_chrom[keep].copy()
    result["_peak_pos"] = on_chrom.index.to_numpy()[keep]
    result["distance_bp"] = distance_bp[keep]
    return result.reset_index(drop=True)


def _find_overlapping_gene_windows(peaks: pd.DataFrame, gene_windows: pd.DataFrame) -> list[list[str]]:
    """For each peak (in row order), which OTHER genes' own promoter-extended
    windows (`core.gene_activity._load_protein_coding_gene_windows` -- the
    same window gene activity uses) it overlaps -- a confound flag, not a
    change to what gets correlated. A "distal" link to the queried gene can
    still just be some OTHER nearby gene's own regulatory activity (a real
    case found on real PBMC data: three of CD14's candidate "links" turned
    out to sit inside CYSTM1/PFDN1's own gene bodies 300kb+ away -- plausibly
    just those genes' own open chromatin, not CD14-specific regulation).
    `gene_windows` should already exclude the queried gene itself (its own
    window was already excluded from `peaks` by `_select_distal_candidate_peaks`,
    so it would never match anyway, but excluding it here too avoids relying
    on that).

    Pure function of plain frames so it's unit-testable on synthetic data,
    same convention as the rest of this module.
    """
    peaks = peaks.reset_index(drop=True)
    overlaps: list[list[str]] = [[] for _ in range(len(peaks))]
    for chrom, gsub in gene_windows.groupby("chrom"):
        psub = peaks[peaks["chrom"] == chrom]
        if psub.empty or gsub.empty:
            continue
        p_start, p_end = psub["start"].to_numpy(), psub["end"].to_numpy()
        p_positions = psub.index.to_numpy()
        g_wstart, g_wend = gsub["win_start"].to_numpy(), gsub["win_end"].to_numpy()
        g_names = gsub["gene_name"].to_numpy()
        overlap = (p_start[:, None] < g_wend[None, :]) & (p_end[:, None] > g_wstart[None, :])
        for i, pos in enumerate(p_positions):
            hits = g_names[overlap[i]]
            if len(hits):
                overlaps[pos] = list(hits)
    return overlaps


def _correlate_candidate_peaks(
    candidates: pd.DataFrame, atac_counts, expr: np.ndarray, gene: str, window_bp: int, top_n: int,
    max_padj: float, min_abs_rho: float,
) -> dict:
    """Spearman-correlate each candidate peak's per-cell accessibility
    (column `_peak_pos` of `atac_counts`) against `expr`, BH-correct across
    all candidates tested for this gene, and keep the top `top_n` by |rho|.

    A peak with zero accessibility variance across these cells (common for
    sparse/low-count distal peaks, especially on a small subsample) makes
    `spearmanr` return `nan` for both rho and pvalue -- `false_discovery_control`
    raises on any NaN input, so such peaks are skipped before correlating
    (undefined correlation, not a real zero), not silently fed in to crash
    the whole gene's lookup. `n_skipped_degenerate_peaks` reports how many,
    so a gene dominated by degenerate peaks is visibly different from one
    with few candidates to begin with.

    "significant" requires BOTH `qvalue < max_padj` AND `|rho| >= min_abs_rho`
    -- q-value alone isn't a safe bar here: at thousands of cells, even a
    trivial rho (~0.1) reaches astronomically small p/q-values (verified on
    the real PBMC fixed-core cache: CD14 had 55/78 candidate peaks clear
    q<0.05 with the strongest rho only 0.235), so a pure significance gate
    would call most of a gene's genomic neighborhood "linked." The effect-
    size floor is what actually distinguishes a real link from large-n noise.

    Each link also carries `overlapping_genes` (from `candidates`, if that
    column is present -- empty list otherwise): OTHER genes whose own
    promoter-extended window this candidate peak sits inside. A "distal"
    link to the queried gene can still just be some other nearby gene's own
    regulatory activity (real example: three of CD14's candidate links on
    the real PBMC data sit inside CYSTM1/PFDN1's own gene bodies 300kb+
    away) -- this is a confound flag only, it never changes which gene's
    expression gets correlated.
    """
    cand_records = candidates.to_dict("records")
    tested, skipped = [], 0
    rhos, pvals = [], []
    for cand, pos in zip(cand_records, candidates["_peak_pos"]):
        col = atac_counts[:, pos]
        col = np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
        if np.ptp(col) == 0:
            skipped += 1
            continue
        rho, pval = spearmanr(col, expr)
        tested.append(cand)
        rhos.append(float(rho))
        pvals.append(float(pval))

    qvals = false_discovery_control(pvals, method="bh") if pvals else []
    links = [
        {
            "chrom": cand["chrom"], "start": int(cand["start"]), "end": int(cand["end"]),
            "distance_bp": int(cand["distance_bp"]), "spearman_rho": rho, "pvalue": pval,
            "qvalue": float(qval), "significant": bool(qval < max_padj and abs(rho) >= min_abs_rho),
            "overlapping_genes": cand.get("overlapping_genes", []),
        }
        for cand, rho, pval, qval in zip(tested, rhos, pvals, qvals)
    ]
    links.sort(key=lambda link: abs(link["spearman_rho"]), reverse=True)
    significant_links = [link for link in links if link["significant"]]
    best_link = max(significant_links, key=lambda link: abs(link["spearman_rho"])) if significant_links else None

    if skipped:
        logger.info("peak_to_gene_links: %s -> skipped %d zero-variance candidate peak(s)", gene, skipped)
    logger.info(
        "peak_to_gene_links: %s -> %d candidate distal peaks, %d significant (best rho=%s)",
        gene, len(links), len(significant_links), best_link["spearman_rho"] if best_link else None,
    )
    return {
        "gene": gene, "window_bp": window_bp, "n_candidate_distal_peaks": len(links),
        "n_skipped_degenerate_peaks": skipped,
        "links": links[:top_n], "best_link": best_link, "any_significant_distal_link": best_link is not None,
    }


def peak_to_gene_links(
    mdata: MuData, gene: str, window_bp: int = DEFAULT_WINDOW_BP,
    proximal_upstream_bp: int = DEFAULT_PROXIMAL_UPSTREAM_BP, top_n: int = 10, max_padj: float = 0.05,
    min_abs_rho: float = DEFAULT_MIN_ABS_RHO, cell_line: str | None = None,
) -> dict:
    """For `gene`, Spearman-correlate each DISTAL candidate peak's per-cell
    accessibility against the gene's own per-cell RNA expression, BH-correct
    across however many candidates were tested, and report the top `top_n`
    by |rho|. Returns enough detail to cite even when nothing links
    (`best_link` is None, `any_significant_distal_link` is False) -- a real
    "not explained, just explainable" result is still a reportable one.

    A link counts as "significant" only if it clears BOTH `max_padj` (BH
    q-value) AND `min_abs_rho` -- q-value alone isn't a safe bar at real
    single-cell sample sizes (see `_correlate_candidate_peaks`'s docstring).

    `cell_line`, if given, restricts BOTH the correlation AND the
    zero-variance check to only that cell line's own cells (matched via
    `rna.obs["cell_line_name"]`) -- a peak invariant WITHIN one line but
    variable pooled across several must correctly read as zero-variance
    there, not inherit variance that only exists between lines. Default
    `None` means pooled across every cell in `mdata`, exactly as before
    this parameter existed -- the right default for a dataset with no
    cell-line concept at all (e.g. the public PBMC data), and the WRONG
    default for a pooled multi-cell-line dataset: a real case caught by
    this project's own adversarial Judger found a genome-wide pooled
    correlation that was actually just two genes both marking the same
    cell line's cluster, not a real within-line relationship (see
    PROGRESS_phase2.md). `n_cells` in the result always reports the real
    cell count actually correlated over, honestly, however small (e.g. a
    164-cell line) -- never silently skipped for being underpowered.

    RNA side uses whatever state `rna.X` is currently in (same caveat as
    `motif_deviations.tf_expression_motif_correlation`: call this after
    clustering's `normalize_total`/`log1p`, not before raw counts). ATAC
    side uses `atac.layers["counts"]` if present (raw accessibility, not
    whatever transform `cluster_atac`'s TF-IDF left in `.X`), falling back
    to `.X` otherwise.
    """
    gene = gene.strip().upper()
    rna = mdata.mod["rna"]
    if gene not in rna.var_names:
        logger.info("peak_to_gene_links: gene %r not found in RNA data", gene)
        return {"error": f"Gene {gene!r} not found in this dataset's RNA var_names."}
    atac = mdata.mod["atac"]

    cell_mask = None
    if cell_line is not None:
        if "cell_line_name" not in rna.obs.columns:
            return {"error": "This dataset has no 'cell_line_name' column -- cell_line scoping isn't applicable here."}
        cell_mask = (rna.obs["cell_line_name"] == cell_line).to_numpy()
        n_cells = int(cell_mask.sum())
        if n_cells == 0:
            return {"error": f"No cells found with cell_line_name == {cell_line!r}."}
    else:
        n_cells = rna.n_obs

    gene_coords = load_protein_coding_gene_coords()
    gene_row = gene_coords[gene_coords["gene_name"] == gene]
    if gene_row.empty:
        logger.info("peak_to_gene_links: gene %r not in protein-coding GENCODE coordinates", gene)
        return {"error": f"Gene {gene!r} has no protein-coding GENCODE coordinates to search peaks around."}
    row = gene_row.iloc[0]

    candidates = _select_distal_candidate_peaks(
        atac.var[["chrom", "start", "end"]], row["chrom"], int(row["start"]), int(row["end"]), row["strand"],
        window_bp=window_bp, proximal_upstream_bp=proximal_upstream_bp,
    )
    if candidates.empty:
        logger.info("peak_to_gene_links: no distal candidate peaks within %d bp of %s", window_bp, gene)
        return {
            "gene": gene, "cell_line": cell_line, "n_cells": n_cells, "window_bp": window_bp,
            "n_candidate_distal_peaks": 0, "n_skipped_degenerate_peaks": 0, "links": [], "best_link": None,
            "any_significant_distal_link": False,
        }

    other_gene_windows = promoter_extend_gene_windows(
        gene_coords[gene_coords["gene_name"] != gene], proximal_upstream_bp,
    )
    candidates = candidates.assign(
        overlapping_genes=_find_overlapping_gene_windows(candidates[["chrom", "start", "end"]], other_gene_windows)
    )

    expr = rna[:, gene].X
    expr = np.asarray(expr.todense()).ravel() if hasattr(expr, "todense") else np.asarray(expr).ravel()
    atac_counts = atac.layers["counts"] if "counts" in atac.layers else atac.X
    if cell_mask is not None:
        expr = expr[cell_mask]
        atac_counts = atac_counts[cell_mask]

    result = _correlate_candidate_peaks(candidates, atac_counts, expr, gene, window_bp, top_n, max_padj, min_abs_rho)
    result["cell_line"] = cell_line
    result["n_cells"] = n_cells
    return result
