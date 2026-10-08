"""Tests for peak-to-gene links (core/peak_to_gene_links.py + the menu
wrapper). The overlap/exclusion and correlation/BH-correction logic is
pure and tested on small synthetic coordinate tables (no real GENCODE
annotation or dataset needed); full end-to-end tests use the shared
session-scoped `agent_fixed_core_mdata` fixture and skip gracefully if raw
data isn't present locally (same convention as test_gene_activity_and_motifs.py).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import spearmanr

from multiome_agent.core.peak_to_gene_links import (
    _correlate_candidate_peaks,
    _find_overlapping_gene_windows,
    _select_distal_candidate_peaks,
    peak_to_gene_links as core_peak_to_gene_links,
)
from multiome_agent.menu.peak_to_gene_links import peak_to_gene_links as menu_peak_to_gene_links


def _six_peaks():
    # Gene: chr1, [100000, 101000), "+". proximal zone (upstream_bp=2000):
    # [98000, 101000). Search window (window_bp=500000): [0, 601000) (clipped
    # at 0 since 100000-500000 is negative).
    return pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr1", "chr2", "chr1"],
        "start": [50000, 99000, 100500, 700000, 50000, 150000],
        "end": [50200, 99500, 100600, 700200, 50200, 150300],
    })  # A=distal-upstream, B=proximal(overlaps prox zone), C=gene-body(in prox zone),
    # D=outside search window, E=wrong chrom, F=distal-downstream


def test_select_distal_candidate_peaks_excludes_proximal_and_out_of_window():
    peaks = _six_peaks()
    result = _select_distal_candidate_peaks(
        peaks, gene_chrom="chr1", gene_start=100000, gene_end=101000, gene_strand="+",
        window_bp=500_000, proximal_upstream_bp=2000,
    )
    # Only A (pos 0) and F (pos 5) survive: B/C fall inside the promoter-extended
    # proximal zone (already covered by gene activity), D is outside the search
    # window, E is on the wrong chromosome.
    assert sorted(result["_peak_pos"].tolist()) == [0, 5]
    by_pos = result.set_index("_peak_pos")
    assert by_pos.loc[0, "distance_bp"] == 100000 - 50200  # upstream of gene start
    assert by_pos.loc[5, "distance_bp"] == 150000 - 101000  # downstream of gene end


def test_select_distal_candidate_peaks_minus_strand_proximal_zone():
    # "-" strand: proximal zone extends downstream of gene end instead of
    # upstream of gene start -- a peak just past gene_end should be excluded,
    # not kept as if it were distal.
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [101000], "end": [101500]})
    result = _select_distal_candidate_peaks(
        peaks, gene_chrom="chr1", gene_start=100000, gene_end=101000, gene_strand="-",
        window_bp=500_000, proximal_upstream_bp=2000,
    )
    assert result.empty


def test_select_distal_candidate_peaks_no_match_on_chrom_returns_empty():
    peaks = pd.DataFrame({"chrom": ["chr2", "chr3"], "start": [1, 1], "end": [100, 100]})
    result = _select_distal_candidate_peaks(
        peaks, gene_chrom="chr1", gene_start=100000, gene_end=101000, gene_strand="+",
    )
    assert result.empty
    assert "_peak_pos" in result.columns and "distance_bp" in result.columns


def test_find_overlapping_gene_windows_flags_peaks_inside_other_genes():
    # peak0 overlaps OTHERGENE's window; peak1 is free (no overlapping gene);
    # peak2 is on a different chromosome than any gene_windows row.
    peaks = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"], "start": [500, 10000, 500], "end": [600, 10100, 600],
    })
    gene_windows = pd.DataFrame({
        "chrom": ["chr1"], "gene_name": ["OTHERGENE"], "win_start": [0], "win_end": [1000],
    })

    overlaps = _find_overlapping_gene_windows(peaks, gene_windows)

    assert overlaps == [["OTHERGENE"], [], []]


def test_find_overlapping_gene_windows_can_flag_multiple_genes_for_one_peak():
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [500], "end": [600]})
    gene_windows = pd.DataFrame({
        "chrom": ["chr1", "chr1"], "gene_name": ["GENEA", "GENEB"],
        "win_start": [0, 400], "win_end": [1000, 700],
    })

    overlaps = _find_overlapping_gene_windows(peaks, gene_windows)

    assert overlaps == [["GENEA", "GENEB"]]


def test_find_overlapping_gene_windows_empty_gene_windows_returns_all_empty():
    peaks = pd.DataFrame({"chrom": ["chr1", "chr2"], "start": [1, 1], "end": [100, 100]})
    gene_windows = pd.DataFrame({"chrom": [], "gene_name": [], "win_start": [], "win_end": []})

    overlaps = _find_overlapping_gene_windows(peaks, gene_windows)

    assert overlaps == [[], []]


def test_correlate_candidate_peaks_ranks_by_abs_rho_and_bh_corrects():
    n_cells = 20
    expr = np.arange(n_cells, dtype=float)
    # peak0: perfectly tracks expression (rho=1). peak1: perfectly anti-tracks
    # (rho=-1, same |rho| as peak0 -- tests that sign doesn't matter for
    # ranking). peak2: no relationship to expression at all.
    atac_counts = np.column_stack([
        np.arange(n_cells, dtype=float),
        -np.arange(n_cells, dtype=float),
        np.array([5, 1, 9, 2, 8, 3, 7, 4, 6, 0, 5, 1, 9, 2, 8, 3, 7, 4, 6, 0], dtype=float),
    ])
    candidates = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1"], "start": [1, 2, 3], "end": [10, 20, 30],
        "distance_bp": [1000, 2000, 3000], "_peak_pos": [0, 1, 2],
    })

    result = _correlate_candidate_peaks(
        candidates, atac_counts, expr, "GENEX", window_bp=500_000, top_n=10, max_padj=0.05, min_abs_rho=0.2,
    )

    assert result["gene"] == "GENEX"
    assert result["n_candidate_distal_peaks"] == 3
    rhos = [link["spearman_rho"] for link in result["links"]]
    assert abs(rhos[0]) == 1.0 and abs(rhos[1]) == 1.0  # the two perfectly-(anti)correlated peaks lead
    assert abs(rhos[2]) < 1.0  # the noise peak trails
    assert result["any_significant_distal_link"] is True
    assert result["best_link"]["significant"] is True
    # `candidates` here has no "overlapping_genes" column at all -- confirms
    # the field still appears on every link, defaulting to an empty list,
    # for callers that haven't computed the confound check.
    assert all(link["overlapping_genes"] == [] for link in result["links"])
    # Independently recompute peak2's rho/pvalue via scipy directly -- the
    # function's own numbers must match, not just be plausible.
    expected_rho, expected_p = spearmanr(atac_counts[:, 2], expr)
    noise_link = next(link for link in result["links"] if link["distance_bp"] == 3000)
    assert noise_link["spearman_rho"] == float(expected_rho)
    assert noise_link["pvalue"] == float(expected_p)


def test_correlate_candidate_peaks_no_significant_link():
    n_cells = 10
    rng_like = np.array([3, 1, 4, 1, 5, 9, 2, 6, 5, 3], dtype=float)
    expr = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1], dtype=float)
    atac_counts = rng_like.reshape(-1, 1)
    candidates = pd.DataFrame({
        "chrom": ["chr1"], "start": [1], "end": [10], "distance_bp": [5000], "_peak_pos": [0],
    })

    result = _correlate_candidate_peaks(
        candidates, atac_counts, expr, "GENEY", window_bp=500_000, top_n=10, max_padj=0.05, min_abs_rho=0.2,
    )

    assert result["any_significant_distal_link"] is False
    assert result["best_link"] is None


def test_correlate_candidate_peaks_skips_zero_variance_peaks_without_crashing():
    # Real bug found during review: spearmanr(constant, anything) returns nan
    # for both rho and pvalue, and scipy's false_discovery_control raises
    # ValueError on any NaN input -- a zero-accessibility distal peak (common
    # for sparse ATAC data) must be skipped before correlating, not crash the
    # whole gene's lookup.
    n_cells = 10
    expr = np.arange(n_cells, dtype=float)
    atac_counts = np.column_stack([
        np.zeros(n_cells),  # degenerate: zero variance -> must be skipped
        np.arange(n_cells, dtype=float),  # real signal
    ])
    candidates = pd.DataFrame({
        "chrom": ["chr1", "chr1"], "start": [1, 2], "end": [10, 20],
        "distance_bp": [1000, 2000], "_peak_pos": [0, 1],
    })

    result = _correlate_candidate_peaks(
        candidates, atac_counts, expr, "GENEZ", window_bp=500_000, top_n=10, max_padj=0.05, min_abs_rho=0.2,
    )

    assert result["n_skipped_degenerate_peaks"] == 1
    assert result["n_candidate_distal_peaks"] == 1
    assert result["links"][0]["distance_bp"] == 2000
    assert result["any_significant_distal_link"] is True


def test_correlate_candidate_peaks_effect_size_floor_overrides_tiny_pvalue():
    # The real motivation for min_abs_rho: at large n, a trivially weak rho
    # still clears q<0.05 easily (exactly what the real CD14 PBMC run showed --
    # see PROGRESS_phase2.md). Build a peak with a real, significant-by-
    # q-value but weak rho (~0.16, below the 0.2 default floor) and confirm
    # it's excluded from "significant" despite a tiny p-value. Fixed seed/
    # scale chosen by direct search so this is a real, reproducible regime,
    # not a guessed one -- rho/p asserted against the live computed values
    # below rather than hardcoded, so the test documents intent, not a magic
    # number that would silently stop testing anything if scipy's RNG stream
    # ever shifted.
    n_cells = 500
    expr = np.arange(n_cells, dtype=float)
    noisy = expr + np.random.RandomState(1).normal(scale=950, size=n_cells)
    candidates = pd.DataFrame({
        "chrom": ["chr1"], "start": [1], "end": [10], "distance_bp": [9000], "_peak_pos": [0],
    })

    result = _correlate_candidate_peaks(
        candidates, noisy.reshape(-1, 1), expr, "GENEW", window_bp=500_000, top_n=10, max_padj=0.05, min_abs_rho=0.2,
    )

    link = result["links"][0]
    assert 0.1 < abs(link["spearman_rho"]) < 0.2  # real, weak-but-nonzero effect
    assert link["qvalue"] < 0.05  # q-value alone would call this "significant"
    assert link["significant"] is False  # the effect-size floor correctly overrides it
    assert result["any_significant_distal_link"] is False


def test_menu_wrapper_unknown_gene_returns_error_not_exception():
    result = menu_peak_to_gene_links(mdata=_FakeMData(), gene="NOT_A_REAL_GENE_XYZ")
    assert "error" in result


def test_peak_to_gene_links_real_data_end_to_end(agent_fixed_core_mdata):
    # CD14 is a real PBMC monocyte marker gene with real GENCODE coordinates and
    # real nearby ATAC peaks -- exercises the full pipeline (real GENCODE parsing,
    # real proximal exclusion, real per-cell correlation + BH correction) against
    # actual data, not just the synthetic unit tests above.
    result = menu_peak_to_gene_links(agent_fixed_core_mdata, "CD14")
    assert "error" not in result
    assert result["gene"] == "CD14"
    assert result["n_candidate_distal_peaks"] >= 0
    assert isinstance(result["any_significant_distal_link"], bool)
    for link in result["links"]:
        assert -1.0 <= link["spearman_rho"] <= 1.0
        assert 0.0 <= link["pvalue"] <= 1.0
        assert 0.0 <= link["qvalue"] <= 1.0
        assert link["distance_bp"] >= 0
    # links must be sorted by |rho| descending
    abs_rhos = [abs(link["spearman_rho"]) for link in result["links"]]
    assert abs_rhos == sorted(abs_rhos, reverse=True)


class _FakeRna:
    var_names = []


class _FakeMData:
    def __init__(self):
        self.mod = {"rna": _FakeRna()}


def _build_cell_line_confound_mdata():
    """2 cell lines (A, B), 10 cells each: GENE1 and a distal peak are both
    simple step functions by line (high in A, low in B) with mutually
    UNCORRELATED within-line jitter -- pooled they look linked (rho=0.714,
    the same numbers verified directly via scipy for the equivalent
    regulon-inference fixture), but within EITHER line alone the
    correlation vanishes (rho=-0.152, n.s.) -- a pure between-line confound,
    not a real distal link.
    """
    rising = np.arange(10, dtype=float)
    noise = np.array([5, 1, 9, 2, 8, 3, 7, 4, 6, 0], dtype=float)
    gene1_expr = np.concatenate([rising + 10, rising])  # line A: 10-19, line B: 0-9
    peak_acc = np.concatenate([noise + 10, noise])  # same between-line step, uncorrelated within-line jitter

    rna = AnnData(X=gene1_expr.reshape(-1, 1), var=pd.DataFrame(index=["GENE1"]))
    rna.obs["cell_line_name"] = ["A"] * 10 + ["B"] * 10

    atac_var = pd.DataFrame({"chrom": ["chr1"], "start": [100000], "end": [100100]}, index=["distal_peak"])
    atac = AnnData(X=peak_acc.reshape(-1, 1), var=atac_var)
    atac.layers["counts"] = peak_acc.reshape(-1, 1)

    class _FakeMuData:
        def __init__(self, rna, atac):
            self.mod = {"rna": rna, "atac": atac}

    return _FakeMuData(rna, atac)


def _gene1_coords():
    return pd.DataFrame({
        "chrom": ["chr1"], "gene_name": ["GENE1"], "strand": ["+"], "start": [10000], "end": [11000],
    })


def test_peak_to_gene_links_pooled_finds_a_spurious_between_line_confound(monkeypatch):
    monkeypatch.setattr("multiome_agent.core.peak_to_gene_links.load_protein_coding_gene_coords", _gene1_coords)
    mdata = _build_cell_line_confound_mdata()

    result = core_peak_to_gene_links(mdata, "GENE1", window_bp=500_000)

    assert result["cell_line"] is None
    assert result["n_cells"] == 20
    assert result["any_significant_distal_link"] is True  # the spurious pooled "signal"


def test_peak_to_gene_links_cell_line_scoping_excludes_the_confound(monkeypatch):
    monkeypatch.setattr("multiome_agent.core.peak_to_gene_links.load_protein_coding_gene_coords", _gene1_coords)
    mdata = _build_cell_line_confound_mdata()

    for line in ("A", "B"):
        result = core_peak_to_gene_links(mdata, "GENE1", window_bp=500_000, cell_line=line)
        assert result["cell_line"] == line
        assert result["n_cells"] == 10  # honest small-n reporting
        assert result["any_significant_distal_link"] is False  # correctly excluded once scoped


def test_peak_to_gene_links_cell_line_no_such_column_returns_error():
    result = menu_peak_to_gene_links(_FakeMData(), "GENE1", cell_line="A")
    assert "error" in result
