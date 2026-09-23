"""Tests for cross-modal cell-type-call validation (core/cross_modal_validation.py).
`best_matching_atac_cluster` is tested on synthetic labels (pure logic, no
real data needed); the full `cross_modal_marker_validation` is tested
against the real, cached PBMC fixed-core result -- known markers should
mostly (not necessarily all -- discordance is expected per CLAUDE.md) cross-
validate.
"""

from __future__ import annotations

import pandas as pd
import pytest
from anndata import AnnData
from mudata import MuData

from multiome_agent.core.cross_modal_validation import best_matching_atac_cluster, cross_modal_marker_validation


def test_best_matching_atac_cluster_picks_majority_overlap():
    # 6 cells: RNA cluster "A" is cells 0-3, all of which land in ATAC
    # cluster "X" (3 of them) or "Y" (1) -- "X" should win by majority overlap.
    rna_labels = pd.Categorical(["A", "A", "A", "A", "B", "B"])
    atac_labels = pd.Categorical(["X", "X", "X", "Y", "Y", "Y"])
    obs_names = [f"cell{i}" for i in range(6)]
    rna = AnnData(obs=pd.DataFrame({"leiden_rna": rna_labels}, index=obs_names))
    atac = AnnData(obs=pd.DataFrame({"leiden_atac": atac_labels}, index=obs_names))
    mdata = MuData({"rna": rna, "atac": atac})

    assert best_matching_atac_cluster(mdata, "A") == "X"
    assert best_matching_atac_cluster(mdata, "B") == "Y"


def test_cross_modal_validation_on_real_pbmc_markers_mostly_confirms(agent_fixed_core_mdata):
    # Uses the cheaper cached agent-tool fixture (n=500, on-disk cache) --
    # not the expensive from-scratch `fixed_core_mdata` (n=3000, rebuilt
    # every session with no cache) -- the validation logic being tested
    # doesn't depend on which real dataset it's checked against.
    results = {
        gene: cross_modal_marker_validation(agent_fixed_core_mdata, gene)
        for gene in ["CD14", "MS4A1", "CD3E", "NKG7"]
    }
    for gene, r in results.items():
        assert r["is_rna_marker"] is True, f"{gene} should be a significant RNA marker of some cluster"
        assert r["matched_atac_cluster"] is not None
        assert isinstance(r["gene_activity_confirms"], bool)
    # Discordance is expected, not a bug (CLAUDE.md) -- at least SOME of
    # these canonical, very-distinct-lineage markers should cross-validate;
    # not requiring all of them (a real, honest mixed result is fine).
    assert any(r["gene_activity_confirms"] for r in results.values())


def test_cross_modal_validation_gene_not_a_marker_returns_none_cleanly(agent_fixed_core_mdata):
    result = cross_modal_marker_validation(agent_fixed_core_mdata, "ZZZNOTAGENE")
    assert result["is_rna_marker"] is False
    assert result["gene_activity_confirms"] is None
