"""Tests for regulon inference (core/regulon_inference.py + the menu
wrapper). The overlap/window-search pure helpers are tested on small
synthetic coordinate tables (no real GENCODE annotation needed); the full
`infer_regulon_targets` orchestration is tested end-to-end on a small,
fully synthetic MuData-like object with `load_protein_coding_gene_coords`
monkeypatched (no real genome annotation or dataset needed for THIS test --
genomic coordinates are fabricated to exercise the promoter-vs-distal
evidence logic precisely); a separate real-data check (not a committed
test) was run manually against the actual PBMC cache once the
motif_match-persisting fixed-core rebuild completed -- see PROGRESS_phase2.md.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from multiome_agent.core.regulon_inference import (
    _genes_near_peaks,
    _overlap_pairs,
    infer_regulon_targets,
)
from multiome_agent.menu.regulon_inference import regulon_inference as menu_regulon_inference


def test_overlap_pairs_finds_peak_gene_window_overlaps():
    peaks = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"], "start": [500, 10000, 500], "end": [600, 10100, 600],
    })
    gene_windows = pd.DataFrame({
        "chrom": ["chr1"], "gene_name": ["GENEA"], "win_start": [0], "win_end": [1000],
    })

    pairs = _overlap_pairs(peaks, gene_windows)

    assert list(pairs["peak_pos"]) == [0]
    assert list(pairs["gene_name"]) == ["GENEA"]


def test_overlap_pairs_empty_when_nothing_overlaps():
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [10000], "end": [10100]})
    gene_windows = pd.DataFrame({"chrom": ["chr1"], "gene_name": ["GENEA"], "win_start": [0], "win_end": [1000]})

    pairs = _overlap_pairs(peaks, gene_windows)

    assert pairs.empty


def test_genes_near_peaks_computes_distance_and_respects_window():
    peaks = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [50000, 500000], "end": [50100, 500100]})
    gene_coords = pd.DataFrame({
        "chrom": ["chr1"], "gene_name": ["GENEB"], "strand": ["+"], "start": [60000], "end": [61000],
    })

    result = _genes_near_peaks(peaks, gene_coords, window_bp=20000)

    # peak0 (distance 9900bp) is within the 20kb window; peak1 (far away) is not.
    assert list(result["peak_pos"]) == [0]
    assert list(result["gene_name"]) == ["GENEB"]
    assert list(result["distance_bp"]) == [9900]


def test_genes_near_peaks_zero_distance_when_overlapping_gene_body():
    peaks = pd.DataFrame({"chrom": ["chr1"], "start": [60500], "end": [60600]})
    gene_coords = pd.DataFrame({
        "chrom": ["chr1"], "gene_name": ["GENEB"], "strand": ["+"], "start": [60000], "end": [61000],
    })

    result = _genes_near_peaks(peaks, gene_coords, window_bp=20000)

    assert list(result["distance_bp"]) == [0]


def _build_synthetic_regulon_mdata():
    """4 candidate genes probing the full promoter-vs-distal evidence logic:
    TARGETA qualifies via direct motif-in-promoter evidence (peak_P);
    TARGETB qualifies via a genuinely distal peak (peak_Q) whose accessibility
    itself significantly tracks TARGETB's expression; TARGETC sits near a
    TF-motif peak (peak_R) too, but that peak's accessibility does NOT track
    TARGETC's expression, so it must be excluded DESPITE TARGETC's own
    TF1-correlation being just as strong as TARGETA/TARGETB's (proving the
    gate filters by real motif+chromatin evidence, not expression alone);
    TARGETD sits near no TF-motif peak at all and is never even considered.
    """
    n_cells = 30
    rising = np.arange(n_cells, dtype=float)
    # Deliberately uncorrelated with `rising` (same pattern used elsewhere in
    # this test suite for a "no real relationship" peak).
    noise = np.array([5, 1, 9, 2, 8, 3, 7, 4, 6, 0] * 3, dtype=float)

    genes = ["TF1", "TARGETA", "TARGETB", "TARGETC", "TARGETD"]
    # All four candidate genes track TF1 equally well in RNA -- the point of
    # this fixture is that only real motif/chromatin evidence (not RNA
    # correlation alone) should admit TARGETC.
    rna_X = np.column_stack([rising] * len(genes))
    rna = AnnData(X=rna_X, var=pd.DataFrame(index=genes))

    peak_names = ["peak_P", "peak_Q", "peak_R"]
    atac_var = pd.DataFrame(
        {"chrom": ["chr1", "chr1", "chr1"], "start": [1000, 50000, 90000], "end": [1100, 50100, 90100]},
        index=peak_names,
    )
    # peak_Q's accessibility tracks `rising` (real distal evidence for
    # TARGETB); peak_R's accessibility is noise (no real distal evidence for
    # TARGETC); peak_P is never correlation-tested (promoter evidence alone
    # qualifies TARGETA).
    atac_X = np.column_stack([noise, rising, noise])
    atac = AnnData(X=atac_X, var=atac_var)
    atac.layers["counts"] = atac_X
    # All 3 peaks carry the TF's motif (column 0 of a single-motif matrix).
    atac.varm["motif_match"] = sp.csr_matrix(np.ones((3, 1), dtype=np.uint8))
    atac.uns["motif_match_names"] = ["MA0001.1.TF1"]

    class _FakeMuData:
        def __init__(self, rna, atac):
            self.mod = {"rna": rna, "atac": atac}

    return _FakeMuData(rna, atac)


def _fake_gene_coords():
    return pd.DataFrame({
        "chrom": ["chr1"] * 5,
        "gene_name": ["TF1", "TARGETA", "TARGETB", "TARGETC", "TARGETD"],
        "strand": ["+"] * 5,
        "start": [500000, 1050, 60000, 100000, 300000],
        "end": [501000, 2050, 61000, 101000, 301000],
    })


def _build_many_targets_mdata(n_genes=6, n_cells=10):
    """All `n_genes` candidate genes qualify via promoter evidence only (one
    dedicated TF-motif peak per gene, far enough apart to never interact) --
    simplest path, no distal complexity needed for this truncation test.
    Each gene's expression is a distinct cyclic shift of TF1's own ranks,
    which reliably gives each gene a DIFFERENT real Spearman correlation
    with TF1 -- the exact ranking isn't hand-derived (cyclic-shift
    correlation isn't simply monotonic in the shift size); the test below
    determines it from an UNCAPPED call's own output instead of assuming it.
    """
    tf_expr = np.arange(n_cells, dtype=float)
    genes = ["TF1"] + [f"TARGET{i}" for i in range(n_genes)]
    rna_cols = [tf_expr] + [np.roll(tf_expr, i + 1).astype(float) for i in range(n_genes)]
    rna = AnnData(X=np.column_stack(rna_cols), var=pd.DataFrame(index=genes))

    peak_names = [f"peak_{i}" for i in range(n_genes)]
    # Gene i's window is [1000 + i*100000, 1100 + i*100000) in the fixture
    # below; its dedicated peak sits well inside that window.
    starts = [1000 + i * 100_000 + 10 for i in range(n_genes)]
    atac_var = pd.DataFrame(
        {"chrom": ["chr1"] * n_genes, "start": starts, "end": [s + 20 for s in starts]}, index=peak_names,
    )
    atac_X = np.ones((n_cells, n_genes), dtype=float)  # accessibility irrelevant -- promoter evidence needs no correlation test
    atac = AnnData(X=atac_X, var=atac_var)
    atac.layers["counts"] = atac_X
    atac.varm["motif_match"] = sp.csr_matrix(np.ones((n_genes, 1), dtype=np.uint8))
    atac.uns["motif_match_names"] = ["MA0001.1.TF1"]

    class _FakeMuData:
        def __init__(self, rna, atac):
            self.mod = {"rna": rna, "atac": atac}

    gene_coords = pd.DataFrame({
        "chrom": ["chr1"] * (n_genes + 1),
        "gene_name": ["TF1"] + [f"TARGET{i}" for i in range(n_genes)],
        "strand": ["+"] * (n_genes + 1),
        "start": [10_000_000] + [1000 + i * 100_000 for i in range(n_genes)],
        "end": [10_001_000] + [1100 + i * 100_000 for i in range(n_genes)],
    })
    return _FakeMuData(rna, atac), gene_coords


def test_infer_regulon_targets_truncates_but_target_gene_always_included(monkeypatch):
    mdata, gene_coords = _build_many_targets_mdata(n_genes=6, n_cells=10)
    monkeypatch.setattr("multiome_agent.core.regulon_inference.load_protein_coding_gene_coords", lambda: gene_coords)

    full = infer_regulon_targets(mdata, "TF1", window_bp=20000, proximal_upstream_bp=2000, top_n=1000)
    assert full["n_candidate_targets"] == 6  # the TRUE count, independent of any cap
    ranked_by_strength = [t["gene"] for t in full["targets"]]  # already sorted strongest-first
    weakest_gene = ranked_by_strength[-1]

    capped = infer_regulon_targets(mdata, "TF1", window_bp=20000, proximal_upstream_bp=2000, top_n=3)
    assert capped["n_candidate_targets"] == 6  # true count unaffected by the cap
    assert len(capped["targets"]) == 3
    assert weakest_gene not in {t["gene"] for t in capped["targets"]}  # truncated away, as expected

    rescued = infer_regulon_targets(
        mdata, "TF1", window_bp=20000, proximal_upstream_bp=2000, top_n=3, target_gene=weakest_gene,
    )
    assert len(rescued["targets"]) == 4  # the 3 strongest + the specifically-requested one
    assert weakest_gene in {t["gene"] for t in rescued["targets"]}
    assert rescued["n_candidate_targets"] == 6  # counts still reflect the true full set


def test_infer_regulon_targets_promoter_and_distal_evidence_both_qualify_excludes_unsupported(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.core.regulon_inference.load_protein_coding_gene_coords", _fake_gene_coords,
    )
    mdata = _build_synthetic_regulon_mdata()

    result = infer_regulon_targets(mdata, "TF1", window_bp=20000, proximal_upstream_bp=2000)

    assert "error" not in result
    assert result["tf"] == "TF1"
    target_genes = {t["gene"] for t in result["targets"]}
    # TARGETC (motif-proximal peak exists, but doesn't itself track TARGETC's
    # expression) and TARGETD (no TF-motif peak nearby at all) must both be
    # excluded, even though TARGETC's RNA-RNA correlation with TF1 would be
    # just as strong as TARGETA/TARGETB's if it had been tested.
    assert target_genes == {"TARGETA", "TARGETB"}

    by_gene = {t["gene"]: t for t in result["targets"]}
    assert by_gene["TARGETA"]["significant"] is True
    assert len(by_gene["TARGETA"]["promoter_motif_evidence"]) == 1
    assert by_gene["TARGETA"]["distal_peak_evidence"] == []

    assert by_gene["TARGETB"]["significant"] is True
    assert by_gene["TARGETB"]["promoter_motif_evidence"] == []
    assert len(by_gene["TARGETB"]["distal_peak_evidence"]) == 1
    assert by_gene["TARGETB"]["distal_peak_evidence"][0]["distance_bp"] == 9900

    assert result["any_significant_target"] is True
    assert result["n_candidate_targets"] == 2


def test_infer_regulon_targets_reports_detection_rate_not_just_marker_status(monkeypatch):
    # Real bug this guards against: a struck-down novelty finding (MYCN->MNX1
    # in SJSA1, see PROGRESS_phase2.md) was originally judged an artifact
    # using cluster-marker status alone ("not a marker of this line's
    # cluster"), when the actually-decisive number -- detection rate within
    # the scoped cell set -- had never been checked. `rising` is 0..29, so
    # exactly 1/30 cells (value 0.0) have zero expression: detection rate
    # must be 29/30 for the TF and for every qualifying target, computed
    # from the SAME (here unscoped) expression arrays already used for the
    # correlation -- not hardcoded or guessed.
    monkeypatch.setattr(
        "multiome_agent.core.regulon_inference.load_protein_coding_gene_coords", _fake_gene_coords,
    )
    mdata = _build_synthetic_regulon_mdata()

    result = infer_regulon_targets(mdata, "TF1", window_bp=20000, proximal_upstream_bp=2000)

    assert "error" not in result
    assert result["tf_detection_rate"] == pytest.approx(29 / 30)
    by_gene = {t["gene"]: t for t in result["targets"]}
    assert by_gene["TARGETA"]["target_detection_rate"] == pytest.approx(29 / 30)
    assert by_gene["TARGETB"]["target_detection_rate"] == pytest.approx(29 / 30)


def test_infer_regulon_targets_tf_not_in_rna_returns_error():
    mdata = _build_synthetic_regulon_mdata()
    result = infer_regulon_targets(mdata, "NOT_A_REAL_TF_XYZ")
    assert "error" in result


def test_infer_regulon_targets_no_motif_match_data_returns_error():
    mdata = _build_synthetic_regulon_mdata()
    del mdata.mod["atac"].uns["motif_match_names"]
    result = infer_regulon_targets(mdata, "TF1")
    assert "error" in result


def test_infer_regulon_targets_no_matching_motif_returns_error(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.core.regulon_inference.load_protein_coding_gene_coords", _fake_gene_coords,
    )
    mdata = _build_synthetic_regulon_mdata()
    result = infer_regulon_targets(mdata, "TARGETA")  # TARGETA is a real RNA gene but has no matching motif
    assert "error" in result


def test_menu_wrapper_unknown_tf_returns_error_not_exception():
    result = menu_regulon_inference(mdata=_build_synthetic_regulon_mdata(), tf_gene="NOT_A_REAL_TF_XYZ")
    assert "error" in result


def _build_cell_line_confound_mdata():
    """2 cell lines (A, B), 10 cells each, demonstrating the EXACT failure mode
    `cell_line` scoping exists to fix: TF1 and CONFOUND are both simple step
    functions by line (high in A, low in B) with mutually UNCORRELATED
    within-line jitter -- so pooled across both lines they look strongly
    correlated (rho=0.714, verified directly via scipy before writing this
    fixture), but within EITHER line alone that correlation vanishes
    (rho=-0.152, n.s.) -- a pure between-line confound, not a real
    relationship. REAL_TARGET, by contrast, tracks TF1's exact rank order
    WITHIN each line too (rho=1.0 both pooled and within either line) -- a
    genuinely real relationship that scoping must NOT destroy.
    """
    rising = np.arange(10, dtype=float)
    noise = np.array([5, 1, 9, 2, 8, 3, 7, 4, 6, 0], dtype=float)
    tf1 = np.concatenate([rising + 10, rising])  # line A: 10-19, line B: 0-9
    confound = np.concatenate([noise + 10, noise])  # same between-line step, uncorrelated within-line jitter
    real_target = np.concatenate([rising + 10, rising])  # tracks TF1 exactly, pooled AND within each line

    genes = ["TF1", "CONFOUND", "REAL_TARGET"]
    rna_X = np.column_stack([tf1, confound, real_target])
    rna = AnnData(X=rna_X, var=pd.DataFrame(index=genes))
    rna.obs["cell_line_name"] = ["A"] * 10 + ["B"] * 10

    peak_names = ["peak_confound_promoter", "peak_real_promoter"]
    atac_var = pd.DataFrame(
        {"chrom": ["chr1", "chr1"], "start": [1000, 50000], "end": [1100, 50100]}, index=peak_names,
    )
    atac_X = np.ones((20, 2), dtype=float)  # accessibility irrelevant here -- promoter evidence needs no correlation
    atac = AnnData(X=atac_X, var=atac_var)
    atac.layers["counts"] = atac_X
    atac.varm["motif_match"] = sp.csr_matrix(np.ones((2, 1), dtype=np.uint8))
    atac.uns["motif_match_names"] = ["MA0001.1.TF1"]

    class _FakeMuData:
        def __init__(self, rna, atac):
            self.mod = {"rna": rna, "atac": atac}

    return _FakeMuData(rna, atac)


def _confound_gene_coords():
    return pd.DataFrame({
        "chrom": ["chr1"] * 3,
        "gene_name": ["TF1", "CONFOUND", "REAL_TARGET"],
        "strand": ["+"] * 3,
        # CONFOUND's promoter-extended window ([1050-2000, 2050] -> clipped [0,2050])
        # contains peak_confound_promoter (1000-1100); REAL_TARGET's window
        # ([50000-2000, 51000] = [48000,51000]) contains peak_real_promoter (50000-50100).
        "start": [500000, 1050, 50000],
        "end": [501000, 2050, 51000],
    })


def test_infer_regulon_targets_pooled_finds_a_spurious_between_line_confound(monkeypatch):
    # Document the EXISTING (pre-fix) failure mode first: pooled across both
    # lines, CONFOUND looks just as "real" as REAL_TARGET.
    monkeypatch.setattr("multiome_agent.core.regulon_inference.load_protein_coding_gene_coords", _confound_gene_coords)
    mdata = _build_cell_line_confound_mdata()

    result = infer_regulon_targets(mdata, "TF1", window_bp=20000, proximal_upstream_bp=2000)

    assert result["cell_line"] is None
    assert result["n_cells"] == 20
    by_gene = {t["gene"]: t for t in result["targets"]}
    assert by_gene["CONFOUND"]["significant"] is True  # the spurious pooled "signal"
    assert by_gene["REAL_TARGET"]["significant"] is True


def test_infer_regulon_targets_cell_line_scoping_excludes_the_confound(monkeypatch):
    monkeypatch.setattr("multiome_agent.core.regulon_inference.load_protein_coding_gene_coords", _confound_gene_coords)
    mdata = _build_cell_line_confound_mdata()

    for line in ("A", "B"):
        result = infer_regulon_targets(mdata, "TF1", window_bp=20000, proximal_upstream_bp=2000, cell_line=line)
        assert result["cell_line"] == line
        assert result["n_cells"] == 10  # honest small-n reporting, not silently inflated or skipped
        by_gene = {t["gene"]: t for t in result["targets"]}
        assert by_gene["CONFOUND"]["significant"] is False  # correctly excluded once scoped
        assert by_gene["REAL_TARGET"]["significant"] is True  # the real relationship survives scoping


def test_infer_regulon_targets_cell_line_no_such_column_returns_error():
    mdata = _build_synthetic_regulon_mdata()  # this fixture has no cell_line_name column at all
    result = infer_regulon_targets(mdata, "TF1", cell_line="A")
    assert "error" in result


def test_infer_regulon_targets_cell_line_no_matching_cells_returns_error():
    mdata = _build_cell_line_confound_mdata()
    result = infer_regulon_targets(mdata, "TF1", cell_line="NOT_A_REAL_LINE")
    assert "error" in result
