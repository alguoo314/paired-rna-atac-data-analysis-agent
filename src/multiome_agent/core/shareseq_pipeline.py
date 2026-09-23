"""Fixed-core-equivalent pipeline for the shareseq-multi-cell-lines multiome
data (see `multiome_agent.data.shareseq_loader`). Deliberately NOT a call
into the tenx-cell-ranger-data `core.pipeline.run_fixed_core` -- that pipeline
assumes raw RNA counts and a real ATAC fragments file, neither of which
this data provides (RNA `.X` is already normalized with no raw layer; ATAC
has a peak-count matrix only, no fragments file), so reusing it as-is would
silently compute wrong numbers rather than fail loudly. Reuses
`core.clustering`'s Leiden/ARI helpers directly where the data shape is
actually compatible (RNA HVG embedding, ATAC LSI, cross-modal ARI).
"""

from __future__ import annotations

import scanpy as sc
from mudata import MuData
from sklearn.metrics import adjusted_rand_score

from multiome_agent.core.clustering import ATAC_CLUSTER_KEY, RNA_CLUSTER_KEY, cluster_atac
from multiome_agent.core.gene_activity import compute_gene_activity_from_peaks
from multiome_agent.core.motif_deviations import compute_motif_deviations
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

# Peaks present in <1% of cells contribute mostly noise to TF-IDF/LSI and
# make it impractically slow at the raw ~1.58M-peak scale (55s LSI at 134K
# peaks after this filter vs. minutes-plus unfiltered, measured directly,
# not assumed) -- a coarser filter than the tenx-cell-ranger pipeline's fixed
# `min_cells=3` (which barely filters anything at this peak count: 1.58M ->
# 1.49M), chosen because this dataset's peak set is a much larger,
# unfiltered/genome-wide set to begin with, not the tenx-cell-ranger data's already
# curated ~100K peaks.
ATAC_MIN_CELL_FRACTION = 0.01


def run_shareseq_qc(mdata: MuData) -> MuData:
    """Surface the real, already-computed QC this data ships with -- no
    recomputation, since `.X` isn't raw counts on the RNA side and there's
    no fragments file on the ATAC side (see module docstring)."""
    rna, atac = mdata.mod["rna"], mdata.mod["atac"]
    logger.info(
        "Shareseq RNA QC (precomputed, reused as-is): median genes/cell=%.0f, "
        "median UMIs/cell=%.0f, median pct_mt=%.2f, doublet rate=%.3f",
        rna.obs["n_genes_by_counts"].median(), rna.obs["total_counts"].median(),
        rna.obs["pct_counts_mt"].median(), (rna.obs["rna_doublet_class"] == "doublet").mean(),
    )
    logger.info(
        "Shareseq ATAC QC (precomputed, reused as-is): median fragments/cell=%.0f, "
        "median FRiP=%.2f, median TSS enrichment=%.1f",
        atac.obs["n_fragment"].median(), atac.obs["frip"].median(), atac.obs["tsse"].median(),
    )
    return mdata


def cluster_shareseq_rna(mdata: MuData, resolution: float = 1.0) -> MuData:
    """Leiden on the precomputed HVG embedding (already normalized + HVG-
    selected upstream -- `core.clustering.cluster_rna` assumes raw counts
    and would double-normalize if called on this data, so this is a
    deliberately separate, shorter pipeline: PCA -> neighbors -> Leiden
    only, no normalize_total/log1p/HVG-selection steps). Cluster labels and
    markers are computed on `rna_hvg` then copied onto the all-genes `rna`
    AnnData so marker lookups can use the full gene set.
    """
    rna_hvg = mdata.mod["rna_hvg"]
    sc.pp.pca(rna_hvg)
    sc.pp.neighbors(rna_hvg)
    sc.tl.leiden(rna_hvg, resolution=resolution, key_added=RNA_CLUSTER_KEY, flavor="igraph")

    rna = mdata.mod["rna"]
    rna.obs[RNA_CLUSTER_KEY] = rna_hvg.obs[RNA_CLUSTER_KEY].reindex(rna.obs_names)
    rna.obs[RNA_CLUSTER_KEY] = rna.obs[RNA_CLUSTER_KEY].astype(rna_hvg.obs[RNA_CLUSTER_KEY].dtype)
    sc.tl.rank_genes_groups(rna, groupby=RNA_CLUSTER_KEY, method="wilcoxon")

    logger.info("Shareseq RNA clustering: %d cells -> %d Leiden clusters", rna.n_obs, rna.obs[RNA_CLUSTER_KEY].nunique())
    return mdata


def cluster_shareseq_atac(mdata: MuData, resolution: float = 1.0) -> MuData:
    """TF-IDF -> LSI -> Leiden on the (peak-filtered) shareseq ATAC matrix,
    reusing `core.clustering.cluster_atac` directly -- ATAC `.X` here IS raw
    counts, so the tenx-cell-ranger pipeline's assumptions hold. Peak-filtered first
    (see `ATAC_MIN_CELL_FRACTION`) since this dataset's peak set is ~15x
    larger than the tenx-cell-ranger data's.
    """
    atac = mdata.mod["atac"]
    sc.pp.filter_genes(atac, min_cells=max(3, int(ATAC_MIN_CELL_FRACTION * atac.n_obs)))
    cluster_atac(atac, resolution=resolution)
    return mdata


def cross_modal_agreement(mdata: MuData) -> dict:
    """RNA-vs-ATAC Leiden cluster agreement (ARI), same metric as the tenx-cell-ranger
    pipeline's `core.clustering.rna_atac_cluster_agreement` -- reimplemented
    rather than called directly since that function expects the labels on
    `mdata["rna"]`/`mdata["atac"]` and here RNA's cluster label technically
    originates from `rna_hvg` before being copied over (harmless either way
    since `rna.obs[RNA_CLUSTER_KEY]` already carries it, but kept explicit)."""
    rna_labels = mdata.mod["rna"].obs[RNA_CLUSTER_KEY]
    atac_labels = mdata.mod["atac"].obs[ATAC_CLUSTER_KEY].reindex(rna_labels.index)
    ari = adjusted_rand_score(rna_labels, atac_labels)
    logger.info("Shareseq data RNA-vs-ATAC cluster agreement: ARI=%.3f", ari)
    return {"ari": ari}


def cell_line_recovery(mdata: MuData) -> dict:
    """Does unsupervised clustering recover the real, genotype-based
    cell-line identity? This dataset's direct ground-truth analog to the
    tenx-cell-ranger pipeline's marker-gene-based "known-biology recovery" --
    real cell-line labels exist here, so recovery can be measured directly
    as ARI against them, separately for RNA and ATAC. Returns only the
    numbers, never a cluster-to-cell-line-name mapping table -- that table
    would contain real cell-line identity strings, which per CLAUDE.md must
    never end up in anything committed (logs, reports, etc.); callers that
    need the mapping for their own uncommitted interactive use can compute
    it themselves from `mdata.mod["rna"].obs["cell_line_name"]` directly.
    """
    rna = mdata.mod["rna"]
    atac = mdata.mod["atac"]
    rna_ari = adjusted_rand_score(rna.obs[RNA_CLUSTER_KEY], rna.obs["cell_line_name"])
    atac_ari = adjusted_rand_score(atac.obs[ATAC_CLUSTER_KEY], atac.obs["cell_line_name"].reindex(atac.obs_names))
    logger.info(
        "Cell-line recovery (ARI vs. true label): RNA=%.3f, ATAC=%.3f", rna_ari, atac_ari
    )
    return {"rna_ari_vs_true_cell_line": rna_ari, "atac_ari_vs_true_cell_line": atac_ari}


def _sanity_check_peak_sequences(mdata: MuData, n_sample: int = 200) -> None:
    """Genome-build sanity check, run BEFORE trusting motif deviations: pull
    real sequence for a sample of this dataset's (peak-filtered) ATAC peaks
    from hg38.2bit and confirm it's real ACGT content in a biologically
    plausible GC range, not degenerate (all-N, which is what you'd see from
    a chrom-naming or genome-build mismatch). Uses the same
    `_add_peak_seq_from_2bit` the real motif-deviation computation uses, not
    a separate reimplementation, so this check exercises the real code path.
    Raises loudly on failure rather than silently proceeding with bad data.
    """
    from multiome_agent.core.motif_deviations import _add_peak_seq_from_2bit

    atac = mdata.mod["atac"]
    sample = atac[:, : min(n_sample, atac.n_vars)].copy()
    sample = _add_peak_seq_from_2bit(sample)
    seqs = sample.uns["peak_seq"]
    all_bases = "".join(seqs)
    n_frac = all_bases.count("N") / max(1, len(all_bases))
    gc_frac = (all_bases.count("G") + all_bases.count("C")) / max(1, len(all_bases) - all_bases.count("N"))
    logger.info(
        "Shareseq ATAC peak-sequence sanity check (n=%d peaks): %.1f%% N bases, %.1f%% GC "
        "(genome-build/chrom-naming mismatch would show as ~100%% N)",
        len(seqs), n_frac * 100, gc_frac * 100,
    )
    if n_frac > 0.5:
        raise ValueError(
            f"{n_frac:.1%} of sampled peak sequence is 'N' -- likely a genome-build or "
            "chrom-naming mismatch between this dataset's ATAC peaks and hg38.2bit. "
            "Not proceeding with motif deviations on this data."
        )
    if not (0.25 <= gc_frac <= 0.70):
        raise ValueError(
            f"Sampled peak GC content ({gc_frac:.1%}) is outside a biologically plausible "
            "range -- likely wrong genome/coordinates. Not proceeding."
        )


def compute_shareseq_motif_deviations(mdata: MuData) -> MuData:
    """chromVAR-style motif deviations for the shareseq-multi-cell-lines data, reusing the
    tenx-cell-ranger pipeline's `compute_motif_deviations` unmodified: it only reads
    `mdata.mod["atac"]` (needs `.var["chrom"/"start"/"end"]`, populated by
    the loader) and `mdata.mod["rna"]` (for TF-expression lookups elsewhere),
    both of which the shareseq MuData already provides in the same shape --
    no adapter needed. Runs a real genome-build sanity check first (see
    `_sanity_check_peak_sequences`) since this dataset's peak set was never
    validated against hg38.2bit before now.
    """
    _sanity_check_peak_sequences(mdata)
    compute_motif_deviations(mdata)
    return mdata


def compute_shareseq_gene_activity(mdata: MuData) -> MuData:
    """Peak-summation-based gene activity for the shareseq-multi-cell-lines data (see
    `gene_activity.compute_gene_activity_from_peaks`) -- the tenx-cell-ranger
    pipeline's fragments-based `compute_gene_activity` can't run here (no
    fragments file for this data, only a peak-count matrix; see
    PROGRESS.md's Stretch-phase entry), but the underlying gene-window
    reference (`snap.genome.hg38.annotation`) is genome-build data, not
    dataset-derived, so it's reused as-is -- no new download needed, only
    the same genome-build sanity check `compute_shareseq_motif_deviations`
    already runs, re-verified here in case this is ever called standalone.
    """
    _sanity_check_peak_sequences(mdata)
    compute_gene_activity_from_peaks(mdata)
    return mdata


def run_shareseq_fixed_core(mdata: MuData) -> MuData:
    """Orchestrates the shareseq-multi-cell-lines pipeline: QC surfacing, clustering
    (both modalities), cross-modal agreement, cell-line recovery, motif
    deviations, gene activity. Stashes results in `.uns`, mirroring the
    tenx-cell-ranger pipeline's convention.
    """
    logger.info("Running shareseq fixed-core-equivalent pipeline on %d cells", mdata.n_obs)
    run_shareseq_qc(mdata)
    cluster_shareseq_rna(mdata)
    cluster_shareseq_atac(mdata)
    mdata.uns["rna_atac_cluster_agreement"] = cross_modal_agreement(mdata)
    mdata.uns["cell_line_recovery"] = cell_line_recovery(mdata)
    compute_shareseq_motif_deviations(mdata)
    compute_shareseq_gene_activity(mdata)
    logger.info("Shareseq fixed-core-equivalent pipeline complete")
    return mdata
