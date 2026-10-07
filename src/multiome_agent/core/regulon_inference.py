"""Regulon inference (TF -> predicted target-gene-set, SCENIC+/RcisTarget-lite):
does a transcription factor's own RNA level track the expression of OTHER,
specific genes it's predicted to regulate -- never an aggregate score, and
never the TF's own expression standing in for a target.

Two-stage, same spirit as real GRN-inference tools (motif + accessibility
define CANDIDATE links; expression correlation is the separate statistic
actually reported):

1. Candidate target genes are genes where the TF's own JASPAR motif is
   found in a nearby accessible peak -- either directly in the gene's own
   promoter-extended window (the standard RcisTarget-style motif-in-cis-
   region criterion), OR in a genuinely distal peak that ALSO clears the
   same significance bar `core.peak_to_gene_links` uses for that specific
   peak/gene pair (so genome-wide proximity to a TF-motif peak alone never
   qualifies a gene -- it has to show a real distal accessibility-expression
   link too).
2. For each candidate (the union of both criteria), the REPORTED statistic
   is always a plain Spearman correlation between the TF's own RNA
   expression and that ONE candidate target gene's own RNA expression --
   deliberately per-target, not one averaged "regulon activity" score,
   because a specific named gene pair (e.g. "SPI1 regulates CD14") is what
   a literature-checking agent can actually look up and confirm, while an
   aggregate module score isn't a citable claim.

Needs `core.motif_deviations.compute_motif_deviations` to have already run
and persisted `atac.varm["motif_match"]` / `atac.uns["motif_match_names"]`
(see that module) -- a dataset's cached fixed-core result built before this
module existed won't have it yet and needs a fresh (not cached) rerun.

Real-data note: some JASPAR motifs are short/degenerate enough to match a
LARGE fraction of all peaks genome-wide (verified on real PBMC data: SPI1's
"MA0080.7.Spi1" ETS-family motif matches 48,557 of 107,353 peaks -- 45%).
For such a motif, testing every genome-wide-nearby (peak, gene) pair one at
a time (a Python-level `scipy.stats.spearmanr` call per pair, plus an O(n)
DataFrame re-scan per hit) made a single TF's regulon take over 10 minutes
without finishing. Fixed by (1) vectorizing the correlation step itself --
rank-transform + z-score each distinct peak/gene row ONCE, then compute
many pairwise Spearman rhos via chunked dot products instead of one scipy
call each (p-values recovered from rho via the same asymptotic t-formula
`scipy.stats.spearmanr` itself uses -- verified numerically identical) --
and (2) capping how many TF-motif peaks the DISTAL search considers for an
over-abundant motif (`DEFAULT_MAX_DISTAL_TF_MOTIF_PEAKS`), reported back as
`distal_search_capped`/`n_tf_motif_peaks_used_for_distal_search` rather than
silently truncated. The PROMOTER-evidence criterion is cheap (a vectorized
overlap test, no per-pair correlation) and stays exhaustive over ALL
TF-motif peaks regardless of this cap.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from mudata import MuData
from scipy.stats import false_discovery_control, rankdata, spearmanr
from scipy.stats import t as t_dist

from multiome_agent.core.gene_activity import load_protein_coding_gene_coords, promoter_extend_gene_windows
from multiome_agent.core.motif_deviations import best_motif_match
from multiome_agent.core.peak_to_gene_links import DEFAULT_MIN_ABS_RHO, DEFAULT_PROXIMAL_UPSTREAM_BP, DEFAULT_WINDOW_BP
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

# How many TF-motif peaks the DISTAL search considers at most -- see module
# docstring's "real-data note." Chosen as a round number comfortably above
# what a genuinely peak-specific motif would ever match, while keeping the
# chunked-correlation step's memory (unique peaks x n_cells, rank-transformed
# float64) in the low hundreds of MB even at this cap.
DEFAULT_MAX_DISTAL_TF_MOTIF_PEAKS = 2000
_CORRELATION_CHUNK_SIZE = 20_000
# How many entries `infer_regulon_targets` returns in `targets` by default --
# see that function's docstring for the real failed-run bug this prevents.
DEFAULT_TARGETS_TOP_N = 30


def _rank_zscore_rows(mat: np.ndarray) -> np.ndarray:
    """Row-wise rank-transform (ties-averaged, matching `scipy.stats.rankdata`'s
    default -- the same ranking `scipy.stats.spearmanr` uses internally) then
    z-score (mean 0, std 1) each row. Pearson correlation between two
    z-scored rank rows equals `spearmanr`'s rho on the original rows, so many
    pairwise Spearman correlations can be computed via plain dot products
    instead of one `scipy.stats.spearmanr` call per pair -- the real fix for
    this module's abundant-motif performance problem (see module docstring).
    Callers must exclude zero-variance rows first (undefined after z-scoring).
    """
    ranks = np.apply_along_axis(rankdata, 1, mat)
    mean = ranks.mean(axis=1, keepdims=True)
    std = ranks.std(axis=1, keepdims=True)
    return (ranks - mean) / std


def _spearman_pvalue_from_rho(rho: np.ndarray, n: int) -> np.ndarray:
    """The same asymptotic t-distribution p-value `scipy.stats.spearmanr`
    itself computes from rho -- verified numerically identical to `spearmanr`'s
    own output on real data before relying on it here. Needed because
    `spearmanr` only returns one pair's p-value per call; batching many pairs
    through `_rank_zscore_rows`'s dot-product trick needs the same formula
    applied vectorized across all of them at once.
    """
    rho = np.clip(rho, -1 + 1e-15, 1 - 1e-15)
    tstat = rho * np.sqrt((n - 2) / (1 - rho**2))
    return 2 * t_dist.sf(np.abs(tstat), df=n - 2)


def _overlap_pairs(peaks: pd.DataFrame, gene_windows: pd.DataFrame) -> pd.DataFrame:
    """Long-form (`peak_pos`, `gene_name`) for every peak/gene-window pair
    that overlaps, grouped by chromosome. `peak_pos` is `peaks`' own row
    POSITION (after an internal `reset_index`), not a pandas index label --
    same convention as the rest of this project's peak/gene overlap helpers
    (see `core.peak_to_gene_links._select_distal_candidate_peaks`). Pure
    function of plain frames, unit-testable on synthetic coordinates.
    """
    peaks = peaks.reset_index(drop=True)
    rows = []
    for chrom, gsub in gene_windows.groupby("chrom"):
        psub = peaks[peaks["chrom"] == chrom]
        if psub.empty or gsub.empty:
            continue
        p_start, p_end = psub["start"].to_numpy(), psub["end"].to_numpy()
        p_positions = psub.index.to_numpy()
        g_wstart, g_wend = gsub["win_start"].to_numpy(), gsub["win_end"].to_numpy()
        g_names = gsub["gene_name"].to_numpy()
        overlap = (p_start[:, None] < g_wend[None, :]) & (p_end[:, None] > g_wstart[None, :])
        pi, gi = np.nonzero(overlap)
        rows.extend(zip(p_positions[pi], g_names[gi]))
    return pd.DataFrame(rows, columns=["peak_pos", "gene_name"])


def _genes_near_peaks(peaks: pd.DataFrame, gene_coords: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    """Long-form (`peak_pos`, `gene_name`, `distance_bp`) for every peak/gene
    pair on the same chromosome within `window_bp` of the gene BODY (0 if
    overlapping) -- the mirror of `core.peak_to_gene_links`'s per-gene
    distal search, but from the peak's side: which genes sit near this one
    TF-motif peak, across the whole genome, rather than which peaks sit
    near one queried gene. `gene_coords` must be
    `load_protein_coding_gene_coords`-shaped (chrom/gene_name/strand/start/end).
    """
    peaks = peaks.reset_index(drop=True)
    rows = []
    for chrom, gsub in gene_coords.groupby("chrom"):
        psub = peaks[peaks["chrom"] == chrom]
        if psub.empty or gsub.empty:
            continue
        p_start, p_end = psub["start"].to_numpy(), psub["end"].to_numpy()
        p_positions = psub.index.to_numpy()
        g_start, g_end = gsub["start"].to_numpy(), gsub["end"].to_numpy()
        g_names = gsub["gene_name"].to_numpy()
        search_start = np.clip(g_start - window_bp, 0, None)
        search_end = g_end + window_bp
        within = (p_start[:, None] < search_end[None, :]) & (p_end[:, None] > search_start[None, :])
        dist = np.where(
            p_start[:, None] > g_end[None, :], p_start[:, None] - g_end[None, :],
            np.where(p_end[:, None] < g_start[None, :], g_start[None, :] - p_end[:, None], 0),
        )
        pi, gi = np.nonzero(within)
        rows.extend(zip(p_positions[pi], g_names[gi], dist[pi, gi].astype(int)))
    return pd.DataFrame(rows, columns=["peak_pos", "gene_name", "distance_bp"])


def infer_regulon_targets(
    mdata: MuData, tf_gene: str, window_bp: int = DEFAULT_WINDOW_BP,
    proximal_upstream_bp: int = DEFAULT_PROXIMAL_UPSTREAM_BP, max_padj: float = 0.05,
    min_abs_rho: float = DEFAULT_MIN_ABS_RHO, target_gene: str | None = None, top_n: int = DEFAULT_TARGETS_TOP_N,
) -> dict:
    """For `tf_gene`, find candidate target genes via real motif + chromatin
    evidence (see module docstring), then Spearman-correlate the TF's own
    RNA against EACH candidate target gene's own RNA (BH-corrected across
    the candidate set). Returns a per-target list -- `targets` -- each with
    its own rho/p/q/significant plus which evidence (promoter and/or
    distal) qualified it, so a downstream literature-checking step gets
    specific, named gene-pair claims rather than one aggregate statistic.

    `targets` is capped to `top_n` entries (all significant ones first, then
    the strongest non-significant ones for context) -- NOT a sample of the
    correlation itself, which is still computed and BH-corrected across
    EVERY real candidate regardless of `top_n`; only the returned LIST is
    bounded. Real bug this fixes: for an abundant motif (e.g. EOMES/RUNX3,
    each matching 15-20% of all peaks), `candidate_genes` can run into the
    thousands, and returning every one of them serialized into a tool
    result blew a real agent conversation past the API's 1M-token prompt
    limit after only a couple of `regulon_inference` calls -- a real failed
    paid run, not a hypothetical. If you already know the specific target
    gene you want to check (e.g. verifying a literature-reported TF-target
    pair), pass `target_gene` -- its full entry is guaranteed to be included
    in `targets` even if it wouldn't otherwise make the `top_n` cut, so
    truncation never hides the one answer you're actually looking for.
    `n_candidate_targets`/`n_significant_targets`/`any_significant_target`
    always reflect the TRUE full counts, never the truncated list's size.

    Returns an error dict (not a raised exception) for a missing gene, a
    missing/unresolvable motif, or missing `motif_match` data (compute_motif_deviations
    must have run and persisted it) -- same convention as `peak_to_gene_links`.
    """
    tf_gene = tf_gene.strip().upper()
    rna = mdata.mod["rna"]
    if tf_gene not in rna.var_names:
        logger.info("infer_regulon_targets: TF %r not found in RNA data", tf_gene)
        return {"error": f"Gene {tf_gene!r} not found in this dataset's RNA var_names."}
    atac = mdata.mod["atac"]

    # `atac.uns["motif_match_names"]` round-trips through h5mu write/read as a
    # numpy array, not the plain Python list it was stored as -- `not <array>`
    # raises ("truth value... is ambiguous") and numpy arrays have no `.index()`
    # method, so normalize to a real list before either is used.
    motif_names = atac.uns.get("motif_match_names")
    motif_names = list(motif_names) if motif_names is not None else None
    if not motif_names or "motif_match" not in atac.varm:
        logger.info("infer_regulon_targets: no motif_match data available")
        return {"error": "No motif_match data available -- compute_motif_deviations must run first."}
    motif_name = best_motif_match(motif_names, tf_gene)
    if motif_name is None:
        logger.info("infer_regulon_targets: no motif matching TF %r", tf_gene)
        return {"error": f"No motif matching {tf_gene!r} among this dataset's motifs."}
    motif_idx = motif_names.index(motif_name)

    # `rna[:, g].X` (AnnData's own single-column indexing) was the dominant
    # cost in a real profiled run (499 of 624 seconds, verified via cProfile
    # on real PBMC data for SPI1's 12,432 candidate targets): single-column
    # slicing out of a CSR matrix is O(total nnz) per call, not O(nnz in that
    # column), and AnnData's view/`_remove_unused_categories` machinery on
    # top of that adds further overhead -- paid thousands of times over for
    # an abundant motif's candidate set. Converting once to CSC (efficient
    # per-column slicing) and indexing the raw matrix directly, bypassing
    # AnnData's `__getitem__` entirely, fixes this at its real source.
    rna_X_csc = rna.X.tocsc() if hasattr(rna.X, "tocsc") else rna.X
    rna_col_idx = {g: i for i, g in enumerate(rna.var_names)}
    rna_cache: dict[str, np.ndarray] = {}

    def gene_expr(g: str) -> np.ndarray:
        if g not in rna_cache:
            col = rna_X_csc[:, rna_col_idx[g]]
            rna_cache[g] = np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
        return rna_cache[g]

    tf_expr = gene_expr(tf_gene)
    if np.ptp(tf_expr) == 0:
        logger.info("infer_regulon_targets: TF %r has zero expression variance", tf_gene)
        return {"error": f"Gene {tf_gene!r} has zero RNA expression variance in this dataset -- cannot correlate."}

    tf_motif_peak_positions = atac.varm["motif_match"][:, motif_idx].nonzero()[0]
    if len(tf_motif_peak_positions) == 0:
        logger.info("infer_regulon_targets: no peaks carry %s's motif %s", tf_gene, motif_name)
        return {
            "tf": tf_gene, "motif": motif_name, "n_tf_motif_peaks": 0, "distal_search_capped": False,
            "n_tf_motif_peaks_used_for_distal_search": 0, "n_candidate_targets": 0,
            "n_significant_targets": 0, "targets": [], "any_significant_target": False,
        }

    tf_peaks_all = atac.var.iloc[tf_motif_peak_positions][["chrom", "start", "end"]].reset_index(drop=True)
    tf_peaks_all["_peak_pos"] = tf_motif_peak_positions
    peak_coord_by_pos = {
        int(pos): (chrom, int(start), int(end))
        for pos, chrom, start, end in zip(
            tf_peaks_all["_peak_pos"], tf_peaks_all["chrom"], tf_peaks_all["start"], tf_peaks_all["end"],
        )
    }

    gene_coords = load_protein_coding_gene_coords()
    other_coords = gene_coords[gene_coords["gene_name"] != tf_gene]
    promoter_windows = promoter_extend_gene_windows(other_coords, proximal_upstream_bp)

    # --- Candidate pass 1: direct motif-in-promoter evidence (no correlation needed to qualify).
    # Cheap (vectorized overlap only) -- stays exhaustive over ALL TF-motif peaks regardless of
    # the distal-search cap below.
    promoter_pairs = _overlap_pairs(tf_peaks_all[["chrom", "start", "end"]], promoter_windows)
    promoter_pairs["atac_pos"] = tf_peaks_all["_peak_pos"].to_numpy()[promoter_pairs["peak_pos"].to_numpy()]
    promoter_evidence: dict[str, list[dict]] = {}
    for atac_pos, gene_name in zip(promoter_pairs["atac_pos"], promoter_pairs["gene_name"]):
        chrom, start, end = peak_coord_by_pos[int(atac_pos)]
        promoter_evidence.setdefault(gene_name, []).append({"chrom": chrom, "start": start, "end": end})

    # --- Candidate pass 2: genuinely distal motif-containing peaks, each requiring its OWN
    # significant peak-accessibility-vs-target-gene-expression link (same bar as peak_to_gene_links)
    # before the gene qualifies -- genome-wide proximity to a TF-motif peak alone is not evidence.
    # Expensive (per-pair correlation), so an over-abundant motif's peak set is capped for THIS
    # search specifically -- see module docstring's "real-data note."
    distal_search_capped = len(tf_peaks_all) > DEFAULT_MAX_DISTAL_TF_MOTIF_PEAKS
    if distal_search_capped:
        sampled_idx = np.random.RandomState(0).choice(
            len(tf_peaks_all), size=DEFAULT_MAX_DISTAL_TF_MOTIF_PEAKS, replace=False,
        )
        tf_peaks_for_distal = tf_peaks_all.iloc[sampled_idx].reset_index(drop=True)
        logger.warning(
            "infer_regulon_targets: %s's motif %s matches %d/%d peaks (too abundant to search "
            "distally in full) -- sampling %d for the distal-evidence search; promoter evidence "
            "stays exhaustive",
            tf_gene, motif_name, len(tf_peaks_all), atac.n_vars, DEFAULT_MAX_DISTAL_TF_MOTIF_PEAKS,
        )
    else:
        tf_peaks_for_distal = tf_peaks_all

    nearby = _genes_near_peaks(tf_peaks_for_distal[["chrom", "start", "end"]], other_coords, window_bp)
    nearby["atac_pos"] = tf_peaks_for_distal["_peak_pos"].to_numpy()[nearby["peak_pos"].to_numpy()]
    if not promoter_pairs.empty:
        nearby = nearby.merge(
            promoter_pairs[["atac_pos", "gene_name"]].assign(_in_promoter=True),
            on=["atac_pos", "gene_name"], how="left",
        )
        nearby = nearby[nearby["_in_promoter"].isna()].drop(columns=["_in_promoter"])
    nearby = nearby[nearby["gene_name"].isin(rna.var_names)]

    # Same CSR-column-slicing cost as `rna_X_csc` above, same fix.
    atac_counts = atac.layers["counts"] if "counts" in atac.layers else atac.X
    atac_counts = atac_counts.tocsc() if hasattr(atac_counts, "tocsc") else atac_counts
    peak_acc_cache: dict[int, np.ndarray] = {}

    def peak_accessibility(pos: int) -> np.ndarray:
        if pos not in peak_acc_cache:
            col = atac_counts[:, pos]
            peak_acc_cache[pos] = np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
        return peak_acc_cache[pos]

    distal_evidence: dict[str, list[dict]] = {}
    if not nearby.empty:
        # Only distinct peaks/genes appearing in `nearby` need rank-transforming, however many
        # (peak, gene) PAIRS reference them -- the real fix for the abundant-motif slowdown.
        unique_pos = [p for p in nearby["atac_pos"].unique() if np.ptp(peak_accessibility(p)) > 0]
        unique_genes = [g for g in nearby["gene_name"].unique() if np.ptp(gene_expr(g)) > 0]
        nearby = nearby[nearby["atac_pos"].isin(unique_pos) & nearby["gene_name"].isin(unique_genes)]
        if not nearby.empty and unique_pos and unique_genes:
            pos_idx = {p: i for i, p in enumerate(unique_pos)}
            gene_idx = {g: i for i, g in enumerate(unique_genes)}
            peak_z = _rank_zscore_rows(np.stack([peak_accessibility(p) for p in unique_pos]))
            gene_z = _rank_zscore_rows(np.stack([gene_expr(g) for g in unique_genes]))
            n_cells = peak_z.shape[1]

            pi = nearby["atac_pos"].map(pos_idx).to_numpy()
            gi = nearby["gene_name"].map(gene_idx).to_numpy()
            rhos = np.empty(len(pi))
            for start in range(0, len(pi), _CORRELATION_CHUNK_SIZE):
                sl = slice(start, start + _CORRELATION_CHUNK_SIZE)
                rhos[sl] = np.einsum("ij,ij->i", peak_z[pi[sl]], gene_z[gi[sl]]) / n_cells
            qvals = false_discovery_control(_spearman_pvalue_from_rho(rhos, n_cells), method="bh")

            sig = np.nonzero((qvals < max_padj) & (np.abs(rhos) >= min_abs_rho))[0]
            nearby_arr = nearby.reset_index(drop=True)
            for i in sig:
                row = nearby_arr.iloc[i]
                chrom, start, end = peak_coord_by_pos[int(row["atac_pos"])]
                distal_evidence.setdefault(row["gene_name"], []).append({
                    "chrom": chrom, "start": start, "end": end, "distance_bp": int(row["distance_bp"]),
                    "spearman_rho": float(rhos[i]), "qvalue": float(qvals[i]),
                })

    candidate_genes = sorted((set(promoter_evidence) | set(distal_evidence)) & set(rna.var_names))
    if not candidate_genes:
        logger.info("infer_regulon_targets: %s -> 0 candidate target genes", tf_gene)
        return {
            "tf": tf_gene, "motif": motif_name, "n_tf_motif_peaks": len(tf_motif_peak_positions),
            "distal_search_capped": distal_search_capped,
            "n_tf_motif_peaks_used_for_distal_search": len(tf_peaks_for_distal),
            "n_candidate_targets": 0, "n_significant_targets": 0, "targets": [], "any_significant_target": False,
        }

    # --- The reported statistic: TF's own RNA vs EACH candidate target gene's own RNA. ---
    tested_genes, trhos, tpvals = [], [], []
    for g in candidate_genes:
        g_expr = gene_expr(g)
        if np.ptp(g_expr) == 0:
            continue
        rho, pval = spearmanr(tf_expr, g_expr)
        tested_genes.append(g)
        trhos.append(float(rho))
        tpvals.append(float(pval))
    tqvals = false_discovery_control(tpvals, method="bh") if tpvals else []

    targets = [
        {
            "gene": g, "spearman_rho": rho, "pvalue": pval, "qvalue": float(qval),
            "significant": bool(qval < max_padj and abs(rho) >= min_abs_rho),
            "promoter_motif_evidence": promoter_evidence.get(g, []),
            "distal_peak_evidence": distal_evidence.get(g, []),
        }
        for g, rho, pval, qval in zip(tested_genes, trhos, tpvals, tqvals)
    ]
    targets.sort(key=lambda t: abs(t["spearman_rho"]), reverse=True)
    n_significant = sum(1 for t in targets if t["significant"])

    # Bound the RETURNED list -- all significant targets first (up to top_n),
    # then the strongest non-significant ones to fill out context; the
    # correlation/BH-correction above still ran over every real candidate
    # regardless. If `target_gene` was asked for specifically and wouldn't
    # otherwise make the cut, append its real entry so truncation never
    # hides the one answer being checked for.
    significant_sorted = [t for t in targets if t["significant"]]
    non_significant_sorted = [t for t in targets if not t["significant"]]
    returned = significant_sorted[:top_n]
    if len(returned) < top_n:
        returned += non_significant_sorted[: top_n - len(returned)]
    if target_gene:
        target_gene_u = target_gene.strip().upper()
        if not any(t["gene"] == target_gene_u for t in returned):
            match = next((t for t in targets if t["gene"] == target_gene_u), None)
            if match:
                returned.append(match)

    logger.info(
        "infer_regulon_targets: %s (%s) -> %d candidate targets, %d significant (%d returned)",
        tf_gene, motif_name, len(targets), n_significant, len(returned),
    )
    return {
        "tf": tf_gene, "motif": motif_name, "n_tf_motif_peaks": len(tf_motif_peak_positions),
        "distal_search_capped": distal_search_capped,
        "n_tf_motif_peaks_used_for_distal_search": len(tf_peaks_for_distal),
        "n_candidate_targets": len(targets), "n_significant_targets": n_significant,
        "targets": returned, "any_significant_target": n_significant > 0,
    }
