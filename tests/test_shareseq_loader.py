"""Tests for the shareseq-multi-cell-lines data loader. Skips gracefully if
`config/local_paths.yaml` isn't set up locally (it's gitignored -- this
data isn't available to everyone who clones this repo).
"""

from __future__ import annotations

import pytest

from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD
from multiome_agent.data.shareseq_loader import load_shareseq_multiome

SHARESEQ_DATA_CONFIGURED = bool(SHARESEQ_RNA_H5AD and SHARESEQ_ATAC_H5AD and SHARESEQ_RNA_HVG_H5AD)

pytestmark = pytest.mark.skipif(
    not SHARESEQ_DATA_CONFIGURED,
    reason="config/local_paths.yaml not set up locally; shareseq-multi-cell-lines data unavailable",
)


@pytest.fixture(scope="module")
def shareseq_mdata():
    return load_shareseq_multiome()


def test_three_modalities_present(shareseq_mdata):
    assert set(shareseq_mdata.mod.keys()) == {"rna", "rna_hvg", "atac"}


def test_barcodes_aligned_across_all_three(shareseq_mdata):
    rna_names = list(shareseq_mdata.mod["rna"].obs_names)
    assert rna_names == list(shareseq_mdata.mod["rna_hvg"].obs_names)
    assert rna_names == list(shareseq_mdata.mod["atac"].obs_names)


def test_rna_hvg_is_gene_subset_of_all_genes(shareseq_mdata):
    hvg_genes = set(shareseq_mdata.mod["rna_hvg"].var_names)
    all_genes = set(shareseq_mdata.mod["rna"].var_names)
    assert hvg_genes.issubset(all_genes)
    assert len(hvg_genes) < len(all_genes)


def test_expected_precomputed_qc_columns_present(shareseq_mdata):
    rna_obs = shareseq_mdata.mod["rna"].obs.columns
    for col in ["total_counts", "pct_counts_mt", "cell_line_name", "rna_doublet_class"]:
        assert col in rna_obs
    atac_obs = shareseq_mdata.mod["atac"].obs.columns
    for col in ["n_fragment", "frip", "tsse", "cell_line_name"]:
        assert col in atac_obs
