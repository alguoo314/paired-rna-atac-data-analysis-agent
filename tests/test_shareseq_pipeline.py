"""Tests for the shareseq-multi-cell-lines fixed-core-equivalent pipeline.
Skips gracefully if the shareseq-multi-cell-lines data isn't configured locally. Asserts
only on numeric outcomes -- never on a literal cell-line/library identity
string, per CLAUDE.md's data-handling rules for this dataset.
"""

from __future__ import annotations

import pytest

from multiome_agent.agent.shareseq_fixed_core_cache import get_shareseq_fixed_core
from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD

SHARESEQ_DATA_CONFIGURED = bool(SHARESEQ_RNA_H5AD and SHARESEQ_ATAC_H5AD and SHARESEQ_RNA_HVG_H5AD)

pytestmark = pytest.mark.skipif(
    not SHARESEQ_DATA_CONFIGURED,
    reason="config/local_paths.yaml not set up locally; shareseq-multi-cell-lines data unavailable",
)


@pytest.fixture(scope="module")
def shareseq_fixed_core_mdata():
    # Uses the disk cache (`agent/shareseq_fixed_core_cache.py`) rather than
    # a fresh `run_shareseq_fixed_core` call -- the pipeline now includes
    # chromVAR motif deviations (~39 min, dominated by background-deviation
    # permutations over ~134K peaks x 5,814 cells x 879 motifs), too slow to
    # recompute per test session. Mirrors the tenx-cell-ranger pipeline's
    # `agent_fixed_core_mdata` cached-fixture pattern in this same conftest.
    return get_shareseq_fixed_core()


def test_rna_clustering_produces_multiple_clusters(shareseq_fixed_core_mdata):
    assert shareseq_fixed_core_mdata.mod["rna"].obs["leiden_rna"].nunique() > 1


def test_atac_clustering_produces_multiple_clusters(shareseq_fixed_core_mdata):
    assert shareseq_fixed_core_mdata.mod["atac"].obs["leiden_atac"].nunique() > 1


def test_cross_modal_agreement_ari_in_valid_range(shareseq_fixed_core_mdata):
    ari = shareseq_fixed_core_mdata.uns["rna_atac_cluster_agreement"]["ari"]
    assert -1.0 <= ari <= 1.0


def test_cell_line_recovery_ari_in_valid_range_and_well_above_chance(shareseq_fixed_core_mdata):
    recovery = shareseq_fixed_core_mdata.uns["cell_line_recovery"]
    assert -1.0 <= recovery["rna_ari_vs_true_cell_line"] <= 1.0
    assert -1.0 <= recovery["atac_ari_vs_true_cell_line"] <= 1.0
    # Real, genotype-confirmed cell-line identity is a strong signal both
    # modalities should recover far above chance-level clustering (ARI~0) --
    # a loose bound (not tuned to the exact observed values, which live only
    # in PROGRESS.md's aggregate numbers) so this doesn't become a brittle
    # regression trap, just a sanity floor.
    assert recovery["rna_ari_vs_true_cell_line"] > 0.3
    assert recovery["atac_ari_vs_true_cell_line"] > 0.3


def test_motif_deviations_matrix_shape_and_sane_values(shareseq_fixed_core_mdata):
    # Correction to an earlier scope decision (see PROGRESS.md): motif
    # deviations only need peak coordinates (parsed from var_names) plus the
    # existing hg38.2bit/JASPAR infra, not the nearest-gene annotation this
    # dataset's peak set genuinely lacks -- those are different things.
    import numpy as np

    dev = shareseq_fixed_core_mdata.mod["atac"].obsm["chromvar_deviations"]
    assert dev.shape[0] == shareseq_fixed_core_mdata.mod["atac"].n_obs
    assert dev.shape[1] > 0
    assert np.issubdtype(np.asarray(dev).dtype, np.floating)  # real numeric matrix, not degenerate/object
    assert not np.isnan(np.asarray(dev)).all()
    assert (dev.std(axis=0) > 0).mean() > 0.5  # most motifs vary across cells, not degenerate/constant


def test_motif_names_recorded(shareseq_fixed_core_mdata):
    names = shareseq_fixed_core_mdata.mod["atac"].uns["chromvar_motif_names"]
    assert len(names) == shareseq_fixed_core_mdata.mod["atac"].obsm["chromvar_deviations"].shape[1]


def test_genome_build_sanity_check_passes_on_real_peaks():
    # Cheap, standalone: no clustering/motif-scan needed, just real peak
    # sequence extraction against hg38.2bit for a small sample -- this is
    # the check that would have caught a genome-build/chrom-naming mismatch
    # before ever starting the expensive full motif-deviation computation.
    from multiome_agent.core.shareseq_pipeline import _sanity_check_peak_sequences
    from multiome_agent.data.shareseq_loader import load_shareseq_multiome

    mdata = load_shareseq_multiome()
    _sanity_check_peak_sequences(mdata, n_sample=50)  # raises on failure


def test_gene_activity_from_peaks_matrix_shape_and_sane_values(shareseq_fixed_core_mdata):
    # Correction to an earlier scope decision: gene activity's fragments-
    # based tenx-cell-ranger method genuinely can't run here (no fragments file), but
    # a peak-summation-based method can, reusing the same GENCODE gene
    # coordinates (genome-build data, not dataset-derived) `compute_gene_activity`
    # already uses -- see `core.gene_activity.compute_gene_activity_from_peaks`.
    import numpy as np

    activity = shareseq_fixed_core_mdata.mod["atac"].obsm["gene_activity"]
    assert activity.shape[0] == shareseq_fixed_core_mdata.mod["atac"].n_obs
    assert activity.shape[1] > 0
    values = activity.toarray() if hasattr(activity, "toarray") else np.asarray(activity)
    assert not np.isnan(values).any()
    assert (values >= 0).all()
    # Real signal, not a degenerate all-zero matrix.
    assert (values.sum(axis=0) > 0).mean() > 0.5


def test_gene_activity_genes_recorded(shareseq_fixed_core_mdata):
    genes = shareseq_fixed_core_mdata.mod["atac"].uns["gene_activity_genes"]
    assert len(genes) == shareseq_fixed_core_mdata.mod["atac"].obsm["gene_activity"].shape[1]


def test_tf_motif_correlation_works_end_to_end_on_shareseq_data():
    # JUN (AP-1 family): broadly expressed, single clean JASPAR match,
    # differentially active across cancer contexts -- the closest analog to
    # SPI1's role in the tenx-cell-ranger pipeline's sanity check, chosen for a
    # cancer-cell-line panel rather than an immune dataset.
    from multiome_agent.core.motif_deviations import tf_expression_motif_correlation

    mdata = get_shareseq_fixed_core()
    result = tf_expression_motif_correlation(mdata, gene="JUN", motif_name_contains="JUN")
    assert "spearman_rho" in result and "pvalue" in result
    assert -1.0 <= result["spearman_rho"] <= 1.0
    assert 0.0 <= result["pvalue"] <= 1.0
