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
    for PBMC (see OWN_DATA_CONTEXT's own reasoning). Per-modality
    signals (each modality's own QC, gene activity, motif deviations, and
    even `cell_line_recovery_ari_atac`) are unaffected by construction --
    the permutation moves each cell's full ATAC row (data AND its own true
    cell_line_name) together, so ATAC stays internally self-consistent;
    only the cross-modal RNA<->ATAC correspondence is actually wrong.
    """
    clean_ari = clean_mdata_with_clusters.uns["rna_atac_cluster_agreement"]["ari"]
    shuffled_ari = cross_modal_agreement(shuffled_mdata)["ari"]
    return {"ari_clean": clean_ari, "ari_faulted": shuffled_ari, "ari_delta": shuffled_ari - clean_ari}


