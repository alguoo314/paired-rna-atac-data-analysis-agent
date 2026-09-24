"""QC-summary builder for the shareseq-multi-cell-lines data, mirroring
`agent/qc_summary.py`'s pattern for the tenx-cell-ranger pipeline. Unlike the
tenx-cell-ranger summary, this one doesn't hardcode one TF's motif correlation (e.g.
SPI1 for PBMC) inline, since there's no single obviously-representative TF
for an 8-cell-line cancer panel the way SPI1 is for PBMC myeloid identity --
`tf_motif_correlation` is available as its own tool for the agent to check
a specific gene on demand (motif deviations ARE computed for this dataset,
see `core/shareseq_pipeline.py`'s `compute_shareseq_motif_deviations`).
Reports aggregate statistics only, per CLAUDE.md's data-handling rules --
no cell-line/library identity strings.
"""

from __future__ import annotations


def shareseq_fixed_core_summary(mdata) -> dict:
    rna, atac = mdata.mod["rna"], mdata.mod["atac"]
    return {
        "n_cells": mdata.n_obs,
        "n_genes": rna.n_vars,
        "n_peaks": atac.n_vars,
        "n_cell_lines": rna.obs["cell_line_name"].nunique(),
        "n_rna_clusters": rna.obs["leiden_rna"].nunique(),
        "n_atac_clusters": atac.obs["leiden_atac"].nunique(),
        "cross_modal_ari": mdata.uns["rna_atac_cluster_agreement"]["ari"],
        "cell_line_recovery_ari_rna": mdata.uns["cell_line_recovery"]["rna_ari_vs_true_cell_line"],
        "cell_line_recovery_ari_atac": mdata.uns["cell_line_recovery"]["atac_ari_vs_true_cell_line"],
        "median_genes_per_cell": float(rna.obs["n_genes_by_counts"].median()),
        "median_umis_per_cell": float(rna.obs["total_counts"].median()),
        "median_pct_mt": float(rna.obs["pct_counts_mt"].median()),
        "median_fragments_per_cell": float(atac.obs["n_fragment"].median()),
        "median_frip": float(atac.obs["frip"].median()),
        "median_tss_enrichment": float(atac.obs["tsse"].median()),
    }


def format_shareseq_qc_summary(summary: dict) -> str:
    return (
        f"{summary['n_cells']} cells pooled from {summary['n_cell_lines']} distinct cell lines, "
        f"{summary['n_genes']} genes, {summary['n_peaks']} ATAC peaks (after feature filtering). "
        f"RNA: {summary['n_rna_clusters']} Leiden clusters, median {summary['median_genes_per_cell']:.0f} "
        f"genes/cell, {summary['median_umis_per_cell']:.0f} UMIs/cell, {summary['median_pct_mt']:.1f}% mito. "
        f"ATAC: {summary['n_atac_clusters']} Leiden clusters, median {summary['median_fragments_per_cell']:.0f} "
        f"fragments/cell, median FRiP {summary['median_frip']:.2f}, median TSS enrichment "
        f"{summary['median_tss_enrichment']:.1f}. "
        f"RNA-ATAC cross-modal cluster agreement (ARI): {summary['cross_modal_ari']:.3f}. "
        f"Recovery of true (genotype-confirmed) cell-line identity from unsupervised clustering "
        f"(ARI vs. ground truth): RNA={summary['cell_line_recovery_ari_rna']:.3f}, "
        f"ATAC={summary['cell_line_recovery_ari_atac']:.3f}. "
        f"chromVAR-style motif accessibility deviations are available for this dataset -- use "
        f"the tf_motif_correlation tool to check a specific transcription factor's expression "
        f"against its own motif's accessibility."
    )
