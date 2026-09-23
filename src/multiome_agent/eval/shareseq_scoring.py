"""Code-level "does this shareseq-multi-cell-lines-data fault actually move the numbers"
signals, mirroring `scoring.py`'s pattern for the tenx-cell-ranger dataset. Reports
statistics only -- no cell-line/library identity strings, per CLAUDE.md's
data-handling rules for this dataset.
"""

from __future__ import annotations

from sklearn.metrics import adjusted_rand_score

from multiome_agent.core.clustering import RNA_CLUSTER_KEY
from multiome_agent.core.shareseq_pipeline import cross_modal_agreement


def score_cell_line_swap(clean_mdata_with_clusters, swapped_mdata) -> dict:
    """Clustering is unsupervised and depends only on RNA/ATAC features, not
    on `cell_line_name` -- so cluster assignments are UNCHANGED by a label
    swap. Cheap to score: reuse the existing clustering, just recompute ARI
    against the corrupted labels rather than rerunning the ~80s pipeline.
    """
    rna = clean_mdata_with_clusters.mod["rna"]
    swapped_labels = swapped_mdata.mod["rna"].obs["cell_line_name"].reindex(rna.obs_names)
    ari_faulted = adjusted_rand_score(rna.obs[RNA_CLUSTER_KEY], swapped_labels)
    ari_clean = clean_mdata_with_clusters.uns["cell_line_recovery"]["rna_ari_vs_true_cell_line"]
    return {"cell_line_recovery_ari_clean": ari_clean, "cell_line_recovery_ari_faulted": ari_faulted}


def score_shuffled_pairing(clean_mdata_with_clusters, shuffled_mdata) -> dict:
    """Cross-modal signal only, mirroring `scoring.score_shuffled_pairing`'s
    tenx-cell-ranger-data version -- but without a TF-motif-rho delta, since there's
    no single canonical TF for an arbitrary cell-line panel the way SPI1 is
    for PBMC (see SHARESEQ_DATASET_CONTEXT's own reasoning). Per-modality
    signals (each modality's own QC, gene activity, motif deviations, and
    even `cell_line_recovery_ari_atac`) are unaffected by construction --
    the permutation moves each cell's full ATAC row (data AND its own true
    cell_line_name) together, so ATAC stays internally self-consistent;
    only the cross-modal RNA<->ATAC correspondence is actually wrong.
    """
    clean_ari = clean_mdata_with_clusters.uns["rna_atac_cluster_agreement"]["ari"]
    shuffled_ari = cross_modal_agreement(shuffled_mdata)["ari"]
    return {"ari_clean": clean_ari, "ari_faulted": shuffled_ari, "ari_delta": shuffled_ari - clean_ari}


def score_mixed_samples(clean_mdata_with_clusters, mixed_mdata, mixed_cell_positions: list[int]) -> dict:
    """Cell-line identity is untouched by this fault (only library/sample
    labels are wrong), so `cell_line_recovery` should NOT move -- a real,
    informative negative result: this fault is invisible to a cell-line-
    identity check by construction, it needs a batch/QC-distribution check
    instead. Scored here via the spread (IQR) of a real per-cell QC metric
    (`pct_counts_mt`) within the now-merged pseudo-sample vs. a real,
    untouched single-library sample of comparable size -- two genuinely
    different technical batches merged into one label should show MORE
    internal QC spread than a real single batch, if the two source
    libraries differ in a QC-relevant property; this is a directional
    expectation to check empirically, not assumed to always hold.
    """
    import numpy as np

    rna = clean_mdata_with_clusters.mod["rna"]
    mixed_rna = mixed_mdata.mod["rna"]

    # The merged pseudo-sample's label is whatever the "kept" library became;
    # every mixed-in cell now carries that same label post-fault.
    merged_label = mixed_rna.obs["library"].iloc[mixed_cell_positions[0]]
    merged_mask = (mixed_rna.obs["library"] == merged_label).to_numpy()
    merged_mt_iqr = float(np.subtract(*np.percentile(mixed_rna.obs["pct_counts_mt"][merged_mask], [75, 25])))

    # Baseline: real single-library samples of comparable size, from the
    # clean (unfaulted) data, for a fair spread comparison.
    clean_lib_sizes = rna.obs["library"].value_counts()
    comparable_libs = clean_lib_sizes[clean_lib_sizes.between(merged_mask.sum() * 0.5, merged_mask.sum() * 1.5)].index
    baseline_iqrs = []
    for lib in comparable_libs:
        vals = rna.obs.loc[rna.obs["library"] == lib, "pct_counts_mt"]
        if len(vals) >= 10:
            baseline_iqrs.append(float(np.subtract(*np.percentile(vals, [75, 25]))))

    return {
        "merged_sample_size": int(merged_mask.sum()),
        "merged_sample_mito_iqr": merged_mt_iqr,
        "comparable_real_single_library_mito_iqr_median": float(np.median(baseline_iqrs)) if baseline_iqrs else None,
        "n_comparable_real_libraries_checked": len(baseline_iqrs),
    }
