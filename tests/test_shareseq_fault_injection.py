"""Tests for the shareseq-multi-cell-lines-only fault types. Skips gracefully if the
shareseq-multi-cell-lines data isn't configured locally. Asserts only on numeric outcomes and
anonymized letters -- never a literal cell-line/library identity string.
"""

from __future__ import annotations

import pytest

from multiome_agent.agent.shareseq_fixed_core_cache import get_shareseq_fixed_core
from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD
from multiome_agent.data.shareseq_loader import load_shareseq_multiome
from multiome_agent.eval.shareseq_fault_injection import shuffle_rna_atac_pairing, swap_cell_line_labels
from multiome_agent.eval.shareseq_scoring import score_cell_line_swap, score_shuffled_pairing

SHARESEQ_DATA_CONFIGURED = bool(SHARESEQ_RNA_H5AD and SHARESEQ_ATAC_H5AD and SHARESEQ_RNA_HVG_H5AD)

pytestmark = pytest.mark.skipif(
    not SHARESEQ_DATA_CONFIGURED,
    reason="config/local_paths.yaml not set up locally; shareseq-multi-cell-lines data unavailable",
)


@pytest.fixture(scope="module")
def clean_mdata():
    return load_shareseq_multiome()


@pytest.fixture(scope="module")
def clean_fixed_core():
    # Cached (see agent/shareseq_fixed_core_cache.py) -- the pipeline now
    # includes chromVAR motif deviations (~39 min), too slow to recompute
    # per test session; these fault-injection tests don't exercise motif
    # deviations at all, but `run_shareseq_fixed_core` always computes them
    # as part of the default pipeline, so avoiding a fresh call here avoids
    # paying that cost for no reason. `.copy()` keeps fault-injection
    # mutations from touching the shared cached object.
    return get_shareseq_fixed_core().copy()


def test_swap_ground_truth_matches_actual_relabeling(clean_mdata):
    faulted, meta = swap_cell_line_labels(clean_mdata, fraction=0.1, seed=0)
    assert meta.fault_type == "cell_line_label_swap"
    n = clean_mdata.n_obs
    assert abs(len(meta.ground_truth["swapped_cell_positions"]) - round(0.1 * n)) <= 1

    true_labels = clean_mdata.mod["rna"].obs["cell_line_name"]
    new_labels = faulted.mod["rna"].obs["cell_line_name"]
    for pos in meta.ground_truth["swapped_cell_positions"]:
        assert new_labels.iloc[pos] != true_labels.iloc[pos]
    # cells NOT in the swap set keep their true label
    unswapped = set(range(n)) - set(meta.ground_truth["swapped_cell_positions"])
    sample_unswapped = list(unswapped)[:50]
    for pos in sample_unswapped:
        assert new_labels.iloc[pos] == true_labels.iloc[pos]
    # ground truth never contains a real cell-line string, only single letters
    for change in meta.ground_truth["label_changes_anonymized"].values():
        assert len(change["true"]) == 1 and change["true"].isalpha()
        assert len(change["fake"]) == 1 and change["fake"].isalpha()


def test_cell_line_swap_degrades_recovery_signal(clean_fixed_core, clean_mdata):
    faulted, _ = swap_cell_line_labels(clean_mdata, fraction=0.3, seed=2)
    result = score_cell_line_swap(clean_fixed_core, faulted)
    assert result["cell_line_recovery_ari_faulted"] < result["cell_line_recovery_ari_clean"]


def test_shuffle_rna_atac_pairing_moves_full_rows_not_just_names(clean_fixed_core):
    faulted, meta = shuffle_rna_atac_pairing(clean_fixed_core, fraction=0.5, seed=1)
    assert meta.fault_type == "shuffled_rna_atac_pairing"
    assert meta.ground_truth["n_shuffled"] > 0
    # rna/rna_hvg untouched; only atac's row order relative to them changed.
    assert list(faulted.mod["rna"].obs_names) == list(clean_fixed_core.mod["rna"].obs_names)
    assert list(faulted.mod["atac"].obs_names) == list(clean_fixed_core.mod["atac"].obs_names)
    # atac's own cell_line_name traveled WITH its permuted row (internally
    # self-consistent), so it should now disagree with rna's for most of the
    # cells whose position was actually shuffled.
    rna_labels = faulted.mod["rna"].obs["cell_line_name"]
    atac_labels = faulted.mod["atac"].obs["cell_line_name"]
    disagreement_rate = (rna_labels.values != atac_labels.values).mean()
    assert disagreement_rate > 0.1  # most (not necessarily all -- a cell could shuffle to a same-cell-line partner)


def test_shuffled_pairing_scoring_degrades_cross_modal_ari(clean_fixed_core):
    faulted, _ = shuffle_rna_atac_pairing(clean_fixed_core, fraction=0.5, seed=1)
    result = score_shuffled_pairing(clean_fixed_core, faulted)
    assert result["ari_faulted"] < result["ari_clean"]
