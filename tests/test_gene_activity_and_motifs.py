"""Tests for fixed-core gene activity scores and chromVAR-style motif
deviations. Uses the shared session-scoped `fixed_core_mdata` fixture from
conftest.py (see its docstring) -- real PBMC data, skips gracefully if raw
data (including hg38.2bit) isn't present locally.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp

from multiome_agent.core.gene_activity import _sum_peaks_into_gene_windows
from multiome_agent.core.motif_deviations import _motif_match_rank


def test_motif_match_rank_prefers_exact_over_substring():
    # Real bug an adversarial "judge" agent caught: "EBF1" substring-matches
    # both the real "MA0154.5.EBF1" and the unrelated "MA0595.1.SREBF1"
    # (contains "ebf1"). Exact token match must rank strictly better.
    assert _motif_match_rank("EBF1", "MA0154.5.EBF1") < _motif_match_rank("EBF1", "MA0595.1.SREBF1")


def test_motif_match_rank_exact_standalone_beats_composite():
    # "CEBPA" is an exact token within the composite "Ddit3::Cebpa", but a
    # real standalone "MA0102.5.CEBPA" motif also exists -- prefer it.
    assert _motif_match_rank("CEBPA", "MA0102.5.CEBPA") < _motif_match_rank("CEBPA", "MA0019.2.Ddit3::Cebpa")


def test_motif_match_rank_composite_beats_pure_substring():
    assert _motif_match_rank("CEBPA", "MA0019.2.Ddit3::Cebpa") < _motif_match_rank("CEBPA", "MA9999.1.NOTCEBPAATALL")


def test_sum_peaks_into_gene_windows_overlap_logic():
    # 4 peaks on chr1, 1 peak on chr2. GENEA's window (0-100) fully contains
    # peak0 and partially overlaps peak1 (which straddles the boundary);
    # peak2 sits entirely outside any window; GENEB's window (80-200)
    # overlaps peak1 too (a deliberately overlapping pair of gene windows,
    # so peak1 should count toward BOTH genes); peak3 is on the wrong
    # chromosome for either gene and must contribute to neither.
    peaks = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr2"],
        "start": [10, 90, 300, 10],
        "end": [50, 110, 350, 50],
    })
    gene_windows = pd.DataFrame({
        "chrom": ["chr1", "chr1"],
        "gene_name": ["GENEA", "GENEB"],
        "win_start": [0, 80],
        "win_end": [100, 200],
    })
    # 2 cells x 4 peaks: cell0 has counts only in peak0/peak1, cell1 only in peak2/peak3.
    X = sp.csr_matrix(np.array([[3, 5, 0, 0], [0, 0, 7, 11]], dtype=np.float64))

    activity, gene_names = _sum_peaks_into_gene_windows(X.tocsc(), peaks, gene_windows)

    assert gene_names == ["GENEA", "GENEB"]
    dense = activity.toarray()
    # GENEA = peak0 + peak1 (peak2/peak3 excluded); GENEB = peak1 only.
    assert dense[0, 0] == 3 + 5  # cell0, GENEA
    assert dense[0, 1] == 5  # cell0, GENEB (peak1 only)
    assert dense[1, 0] == 0  # cell1 has no counts in peak0/peak1
    assert dense[1, 1] == 0


def test_sum_peaks_into_gene_windows_no_overlap_returns_empty():
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [1000], "end": [1050]})
    gene_windows = pd.DataFrame(
        {"chrom": ["chr1"], "gene_name": ["GENEA"], "win_start": [0], "win_end": [100]}
    )
    X = sp.csr_matrix(np.array([[5]], dtype=np.float64))

    activity, gene_names = _sum_peaks_into_gene_windows(X.tocsc(), peaks, gene_windows)

    assert gene_names == []
    assert activity.shape == (1, 0)


def test_gene_activity_matrix_shape_and_sane_values(fixed_core_mdata):
    activity = fixed_core_mdata["atac"].obsm["gene_activity"]
    assert activity.shape[0] == fixed_core_mdata.n_obs
    assert activity.shape[1] > 0
    # `gene_activity` is a scipy sparse matrix (see gene_activity.py docstring
    # for why it's not a pandas sparse-dtype DataFrame), so densify via
    # `.toarray()` rather than pandas' `.sparse` accessor.
    values = np.asarray(activity.toarray() if hasattr(activity, "toarray") else activity)
    assert not np.isnan(values).any()
    assert (values >= 0).all()


def test_gene_activity_genes_recorded(fixed_core_mdata):
    genes = fixed_core_mdata["atac"].uns["gene_activity_genes"]
    assert len(genes) == fixed_core_mdata["atac"].obsm["gene_activity"].shape[1]


def test_motif_deviations_matrix_shape_and_sane_values(fixed_core_mdata):
    dev = fixed_core_mdata["atac"].obsm["chromvar_deviations"]
    assert dev.shape[0] == fixed_core_mdata.n_obs
    assert dev.shape[1] > 0
    assert not dev.isna().all(axis=None)
    # z-scores should vary across cells for at least most motifs, not be
    # uniformly zero (which would indicate the bias-correction step silently
    # produced degenerate output rather than genuinely failing).
    assert (dev.std(axis=0) > 0).mean() > 0.5


def test_motif_names_recorded(fixed_core_mdata):
    names = fixed_core_mdata["atac"].uns["chromvar_motif_names"]
    assert len(names) == fixed_core_mdata["atac"].obsm["chromvar_deviations"].shape[1]
    assert any("spi1" in n.lower() for n in names)


def test_tf_motif_sanity_check_recorded_and_positive(fixed_core_mdata):
    check = fixed_core_mdata.uns["tf_motif_sanity_check"]
    assert check["gene"] == "SPI1"
    # A moderate positive correlation is the expected result (RNA-ATAC
    # relationships are noisy per CLAUDE.md); this is a low bar deliberately
    # -- a sanity check, not a strict pass/fail gate.
    assert check["spearman_rho"] > 0.1
    assert check["pvalue"] < 0.05
