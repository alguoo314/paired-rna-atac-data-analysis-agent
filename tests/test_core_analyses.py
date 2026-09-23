"""Tests for the fixed-core RNA/ATAC QC, clustering, and cluster agreement.

Uses the session-scoped `fixed_core_mdata` fixture from conftest.py (real
PBMC data, run once for the whole test session) -- TSS enrichment alone
takes a few minutes (snapatac2 sorts the whole fragments file per call), so
this is deliberately not re-run per test or per test file. Skips gracefully
if raw data isn't present locally.
"""

import pytest

from multiome_agent.core.clustering import ATAC_CLUSTER_KEY, RNA_CLUSTER_KEY, marker_is_recovered

CANONICAL_PBMC_MARKERS = ["CD14", "LYZ", "MS4A1", "CD3E", "NKG7"]


def test_rna_qc_columns_in_sane_ranges(fixed_core_mdata):
    obs = fixed_core_mdata["rna"].obs
    assert (obs["total_counts"] > 0).all()
    assert (obs["n_genes_by_counts"] > 0).all()
    assert obs["pct_counts_mt"].between(0, 100).all()
    assert obs["doublet_score"].between(0, 1).all()
    assert obs["predicted_doublet"].dtype == bool


def test_atac_qc_columns_in_sane_ranges(fixed_core_mdata):
    obs = fixed_core_mdata["atac"].obs
    assert obs["frip"].between(0, 1).all()
    assert (obs["tss_enrichment"] > 0).all()
    assert (obs["nucleosome_signal"] > 0).all()


def test_rna_clustering_produces_multiple_clusters(fixed_core_mdata):
    assert fixed_core_mdata["rna"].obs[RNA_CLUSTER_KEY].nunique() > 1


def test_atac_clustering_produces_multiple_clusters(fixed_core_mdata):
    assert fixed_core_mdata["atac"].obs[ATAC_CLUSTER_KEY].nunique() > 1


def test_cluster_agreement_ari_in_valid_range(fixed_core_mdata):
    ari = fixed_core_mdata.uns["rna_atac_cluster_agreement"]["ari"]
    assert -1.0 <= ari <= 1.0


def test_cluster_agreement_contingency_table_shape(fixed_core_mdata):
    contingency = fixed_core_mdata.uns["rna_atac_cluster_agreement"]["contingency"]
    rna_n_clusters = fixed_core_mdata["rna"].obs[RNA_CLUSTER_KEY].nunique()
    atac_n_clusters = fixed_core_mdata["atac"].obs[ATAC_CLUSTER_KEY].nunique()
    assert contingency.shape == (rna_n_clusters, atac_n_clusters)
    assert contingency.to_numpy().sum() == fixed_core_mdata.n_obs


@pytest.mark.parametrize("gene", CANONICAL_PBMC_MARKERS)
def test_known_pbmc_marker_recovered(fixed_core_mdata, gene):
    assert marker_is_recovered(fixed_core_mdata["rna"], gene), (
        f"{gene} should be a significant (padj<0.05, logFC>1) marker of at least one cluster"
    )
