"""Shared fixed-core summary builder, used by both the report writer's
"dataset & fixed-core summary" section and the agent's `get_qc_summary` tool
-- one implementation, not two independently-drifting ones.
"""

from __future__ import annotations

from multiome_agent.core.clustering import ATAC_CLUSTER_KEY, RNA_CLUSTER_KEY
from multiome_agent.core.motif_deviations import tf_expression_motif_correlation


def fixed_core_summary(mdata) -> dict:
    rna, atac = mdata.mod["rna"], mdata.mod["atac"]
    ari = mdata.uns["rna_atac_cluster_agreement"]["ari"]
    spi1 = tf_expression_motif_correlation(mdata, gene="SPI1", motif_name_contains="Spi1")
    return {
        "n_cells": mdata.n_obs,
        "n_genes": rna.n_vars,
        "n_peaks": atac.n_vars,
        "n_rna_clusters": rna.obs[RNA_CLUSTER_KEY].nunique(),
        "n_atac_clusters": atac.obs[ATAC_CLUSTER_KEY].nunique(),
        "ari": ari,
        "median_genes_per_cell": float(rna.obs["n_genes_by_counts"].median()),
        "median_umis_per_cell": float(rna.obs["total_counts"].median()),
        "median_pct_mt": float(rna.obs["pct_counts_mt"].median()),
        "predicted_doublet_rate": float(rna.obs["predicted_doublet"].mean()),
        "median_doublet_score": float(rna.obs["doublet_score"].median()),
        "median_frip": float(atac.obs["frip"].median()),
        "median_tss_enrichment": float(atac.obs["tss_enrichment"].median()),
        "median_nucleosome_signal": float(atac.obs["nucleosome_signal"].median()),
        "median_atac_fragments": float(atac.obs["atac_fragments"].median()),
        "spi1_motif_rho": float(spi1["spearman_rho"]),
    }


def format_qc_summary(summary: dict) -> str:
    """Plain-text rendering for the agent's `get_qc_summary` tool result --
    every number a claim in an eval answer should be able to trace back to."""
    return (
        f"{summary['n_cells']} cells, {summary['n_genes']} genes, {summary['n_peaks']} ATAC peaks. "
        f"RNA: {summary['n_rna_clusters']} Leiden clusters, median {summary['median_genes_per_cell']:.0f} "
        f"genes/cell, {summary['median_umis_per_cell']:.0f} UMIs/cell, {summary['median_pct_mt']:.1f}% mito, "
        f"{summary['predicted_doublet_rate']:.1%} predicted doublets (median doublet_score "
        f"{summary['median_doublet_score']:.3f}). "
        f"ATAC: {summary['n_atac_clusters']} Leiden clusters, median {summary['median_atac_fragments']:.0f} "
        f"fragments/cell, median FRiP {summary['median_frip']:.2f}, "
        f"median TSS enrichment {summary['median_tss_enrichment']:.1f}, "
        f"median nucleosome signal {summary['median_nucleosome_signal']:.2f}. "
        f"RNA-ATAC cluster agreement (ARI): {summary['ari']:.3f}. "
        f"SPI1 expression vs. its own motif's chromVAR deviation (Spearman rho): "
        f"{summary['spi1_motif_rho']:.3f}."
    )
