"""Tests for the known-biology checklist scorer. Reuses the cached
`agent_fixed_core_mdata` fixture (real data, real numbers) -- see
conftest.py.
"""

from __future__ import annotations

import pytest

from multiome_agent.eval.checklist_scoring import annotate_clusters, load_checklist, score_checklist


def test_checklist_loads_with_expected_sections():
    checklist = load_checklist()
    assert len(checklist["rna_markers"]) == 5
    assert len(checklist["motifs"]) == 6
    assert len(checklist["tf_expression_tracks_motif"]) == 1


def test_cluster_annotation_produces_expected_cell_types(agent_fixed_core_mdata):
    checklist = load_checklist()
    labels = annotate_clusters(agent_fixed_core_mdata["rna"], checklist)
    labeled_types = set(v for v in labels.values() if v is not None)
    # Every checklist cell type should label at least one cluster in a
    # dataset this size -- if none did, the annotation logic or the
    # underlying clustering would be broken, not just small-sample noise.
    assert {"monocyte", "B cell", "T cell", "NK cell"}.issubset(labeled_types)


@pytest.fixture(scope="module")
def scorecard(agent_fixed_core_mdata):
    return score_checklist(agent_fixed_core_mdata)


def test_all_rna_markers_recovered(scorecard):
    assert scorecard["summary"]["rna_markers"]["recall"] == 1.0
    for r in scorecard["rna_markers"]:
        assert r["recovered"], f"{r['gene']} ({r['cell_type']}) not recovered as a significant marker"


def test_tf_expression_tracks_motif_recovered(scorecard):
    assert scorecard["summary"]["tf_expression_tracks_motif"]["recall"] == 1.0


def test_motif_scorecard_structure_and_reasonable_recall(scorecard):
    # Real biology is noisier at n=500 -- expect most but not necessarily
    # all motifs recovered (PAX5 is a documented, genuine miss at this
    # subsample size, not a scorer bug -- see PROGRESS.md). A >=0.5 floor
    # catches a real regression without being brittle to which specific
    # motif is weak on a given seed.
    assert scorecard["summary"]["motifs"]["recall"] >= 0.5
    for motif in scorecard["motifs"]:
        assert motif["status"] in ("recovered", "not_recovered", "present_but_not_scoreable")
        assert motif["matched_columns"], f"{motif['motif_name']} should exist in the 879-motif JASPAR set"
