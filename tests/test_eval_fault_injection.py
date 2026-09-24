"""Tests for the fault injector and its code-level scoring signals.

Cheap unit tests (ground truth correctness, shapes) run unconditionally.
Tests that check a fault actually shifts fixed-core numbers reuse the
session-scoped `agent_fixed_core_mdata` fixture (n=500, already computed and
cached from steps 1-5 -- see conftest.py) rather than paying for a fresh
~17-minute fixed-core run, and budget real snapatac2 reruns carefully
(downsampling needs one; shuffled-pairing and doublets don't).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from multiome_agent.core.rna_qc import run_rna_qc
from multiome_agent.data.loader import RAW_DIR, load_pbmc_multiome
from multiome_agent.eval.fault_injection import (
    clean_control,
    downsample_atac_fragments,
    inject_doublets,
    shuffle_rna_atac_pairing,
)
from multiome_agent.eval.scoring import score_doublets, score_downsampling, score_shuffled_pairing

FRAGMENTS_FILE = RAW_DIR / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
RAW_DATA_PRESENT = FRAGMENTS_FILE.exists()

pytestmark = pytest.mark.skipif(not RAW_DATA_PRESENT, reason="raw PBMC ATAC fragments not present in data/raw/")


# ---- cheap unit tests: no fixed-core run needed ----


def test_clean_control_is_unmodified():
    mdata = load_pbmc_multiome(n_cells=100, seed=0, force_reload=False)
    cc, meta = clean_control(mdata)
    assert meta.fault_type == "clean_control"
    assert meta.severity == 0.0
    assert cc.n_obs == mdata.n_obs
    assert (cc["rna"].X != mdata["rna"].X).nnz == 0


@pytest.mark.parametrize("fraction", [0.0, 0.2, 0.5, 1.0])
def test_shuffle_ground_truth_matches_actual_permutation(fraction):
    mdata = load_pbmc_multiome(n_cells=200, seed=0, force_reload=False)
    shuffled, meta = shuffle_rna_atac_pairing(mdata, fraction=fraction, seed=1)

    assert meta.severity == fraction
    assert list(shuffled["rna"].obs_names) == list(shuffled["atac"].obs_names) == list(mdata["rna"].obs_names)

    true_source = meta.ground_truth["true_atac_source_by_rna_barcode"]
    for barcode, source_barcode in true_source.items():
        assert barcode != source_barcode
        expected_row = mdata["atac"][source_barcode].X.toarray()
        actual_row = shuffled["atac"][barcode].X.toarray()
        assert np.array_equal(expected_row, actual_row)
    if fraction > 0:
        assert len(true_source) > 0
    if fraction == 1.0:
        # every cell should have moved (no accidental fixed points survive a rotation)
        assert len(true_source) == mdata.n_obs


def test_shuffle_zero_fraction_is_a_true_no_op():
    mdata = load_pbmc_multiome(n_cells=100, seed=0, force_reload=False)
    shuffled, meta = shuffle_rna_atac_pairing(mdata, fraction=0.0, seed=1)
    assert meta.ground_truth["true_atac_source_by_rna_barcode"] == {}
    assert (shuffled["atac"].X != mdata["atac"].X).nnz == 0


def test_inject_doublets_ground_truth_and_count_summation():
    mdata = load_pbmc_multiome(n_cells=200, seed=0, force_reload=False)
    faulted, meta = inject_doublets(mdata, doublet_rate=0.1, seed=2)

    names = meta.ground_truth["doublet_obs_names"]
    pairs = meta.ground_truth["source_pairs"]
    assert len(names) == len(pairs)
    assert faulted.n_obs == mdata.n_obs + len(names)
    assert abs(len(names) / faulted.n_obs - 0.1) < 0.02

    a, b = pairs[0]
    expected = mdata["rna"].X[a].toarray() + mdata["rna"].X[b].toarray()
    actual = faulted["rna"][names[0]].X.toarray()
    assert np.array_equal(expected, actual)
    # doublets aren't real droplets -- no fragments-derived ATAC QC for them
    assert faulted["atac"].obs.loc[names[0], "atac_fragments"] != faulted["atac"].obs.loc[names[0], "atac_fragments"]  # NaN


def test_downsample_fragments_file_keep_frac_is_approximately_respected(tmp_path):
    mdata = load_pbmc_multiome(n_cells=50, seed=0, force_reload=False)
    barcodes = list(mdata["atac"].obs_names)

    full_path, full_meta = downsample_atac_fragments(FRAGMENTS_FILE, barcodes, keep_frac=1.0, out_dir=tmp_path, seed=0)
    half_path, half_meta = downsample_atac_fragments(FRAGMENTS_FILE, barcodes, keep_frac=0.5, out_dir=tmp_path, seed=0)

    assert full_meta.severity == 0.0
    assert half_meta.severity == 0.5
    kept_ratio = half_meta.ground_truth["n_fragments_kept"] / full_meta.ground_truth["n_fragments_kept"]
    assert 0.4 < kept_ratio < 0.6


# ---- signal-shift tests: reuse cached fixed-core results where possible ----


def test_shuffled_pairing_degrades_cross_modal_signals_monotonically(agent_fixed_core_mdata):
    prev_ari = None
    prev_rho = None
    for fraction in [0.2, 0.5, 1.0]:
        shuffled, _ = shuffle_rna_atac_pairing(agent_fixed_core_mdata, fraction=fraction, seed=1)
        result = score_shuffled_pairing(agent_fixed_core_mdata, shuffled)
        assert result["ari_faulted"] <= result["ari_clean"]
        assert abs(result["tf_motif_rho_faulted"]) <= abs(result["tf_motif_rho_clean"])
        if prev_ari is not None:
            assert result["ari_faulted"] <= prev_ari + 0.05  # allow small noise, expect monotonic-ish decline
        prev_ari, prev_rho = result["ari_faulted"], result["tf_motif_rho_faulted"]
    # at 100% shuffle (a near-total derangement), agreement should collapse close to chance
    assert result["ari_faulted"] < 0.05


def test_clean_vs_clean_shuffle_shows_no_spurious_signal(agent_fixed_core_mdata):
    """False-alarm sanity check on the scoring machinery itself: comparing
    the clean result against itself (fraction=0.0) must show exactly zero
    delta, not just a small one."""
    unshuffled, _ = shuffle_rna_atac_pairing(agent_fixed_core_mdata, fraction=0.0, seed=1)
    result = score_shuffled_pairing(agent_fixed_core_mdata, unshuffled)
    assert result["ari_delta"] == pytest.approx(0.0, abs=1e-9)
    assert result["tf_motif_rho_delta"] == pytest.approx(0.0, abs=1e-9)


def test_doublets_elevate_detection_rate_among_synthetic_cells():
    mdata = load_pbmc_multiome(n_cells=500, seed=0, force_reload=False)
    faulted, meta = inject_doublets(mdata, doublet_rate=0.1, seed=3)
    rna = faulted["rna"].copy()
    run_rna_qc(rna)
    result = score_doublets(rna.obs, meta.ground_truth["doublet_obs_names"])
    assert result["doublet_score_median_synthetic"] > result["doublet_score_median_real"]
    assert result["predicted_doublet_rate_synthetic"] > result["predicted_doublet_rate_real"]


@pytest.mark.slow
def test_downsampling_reduces_fragments_derived_atac_signals(tmp_path, agent_fixed_core_mdata):
    barcodes = list(agent_fixed_core_mdata["atac"].obs_names)
    downsampled_path, meta = downsample_atac_fragments(FRAGMENTS_FILE, barcodes, keep_frac=0.2, out_dir=tmp_path, seed=0)

    result = score_downsampling(
        agent_fixed_core_mdata["atac"], barcodes, FRAGMENTS_FILE, downsampled_path, keep_frac=0.2
    )
    print("\ndownsampling (keep_frac=0.2) scoring result:", result)
    assert result["n_fragment_median_faulted"] < result["n_fragment_median_clean"]
    assert result["n_fragment_median_pct_change"] < -50
    assert result["tss_enrichment_median_faulted"] < result["tss_enrichment_median_clean"]
    # nucleosome_signal is a ratio (mono-nucleosomal / nucleosome-free
    # fragments), and `_stream_downsample_fragments` drops lines uniformly at
    # random regardless of fragment length -- so unlike n_fragment/
    # tss_enrichment (which are absolute-depth-sensitive and should clearly
    # drop), there's no length bias to justify asserting a specific
    # direction here. Just check it's still a valid, finite, non-negative
    # number after downsampling -- a sanity check on the `n=` fix, not a
    # directional claim.
    assert np.isfinite(result["nucleosome_signal_median_faulted"])
    assert result["nucleosome_signal_median_faulted"] >= 0
