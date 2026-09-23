"""Cheap tests for the shareseq-multi-cell-lines-data novelty-check splitting logic (no API
cost -- the real 4-agent-run investigation is validated separately and
documented with real numbers in PROGRESS.md, not re-run on every test pass).
"""

from __future__ import annotations

import pytest

from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD
from multiome_agent.eval.shareseq_novelty_checks import split_cells_in_half

SHARESEQ_DATA_CONFIGURED = bool(SHARESEQ_RNA_H5AD and SHARESEQ_ATAC_H5AD and SHARESEQ_RNA_HVG_H5AD)
pytestmark = pytest.mark.skipif(
    not SHARESEQ_DATA_CONFIGURED,
    reason="config/local_paths.yaml not set up locally; shareseq-multi-cell-lines data unavailable",
)


def test_split_is_disjoint_and_covers_all_cells():
    from multiome_agent.data.shareseq_loader import load_shareseq_multiome

    mdata = load_shareseq_multiome()
    half1, half2 = split_cells_in_half(mdata, seed=0)

    assert set(half1).isdisjoint(set(half2))
    assert len(half1) + len(half2) in (mdata.n_obs, mdata.n_obs - 1)  # odd n_obs drops one cell to a clean 50/50
    assert abs(len(half1) - len(half2)) <= 1


def test_split_is_reproducible_with_same_seed():
    from multiome_agent.data.shareseq_loader import load_shareseq_multiome

    mdata = load_shareseq_multiome()
    a1, a2 = split_cells_in_half(mdata, seed=0)
    b1, b2 = split_cells_in_half(mdata, seed=0)
    assert a1 == b1
    assert a2 == b2


def test_different_seeds_give_different_splits():
    from multiome_agent.data.shareseq_loader import load_shareseq_multiome

    mdata = load_shareseq_multiome()
    a1, _ = split_cells_in_half(mdata, seed=0)
    b1, _ = split_cells_in_half(mdata, seed=1)
    assert a1 != b1
