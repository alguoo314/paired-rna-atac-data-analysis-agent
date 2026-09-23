"""Tests for the known-biology checklist generator (agent/checklist_generator.py).

Unit-tests the extraction logic (record_checklist_item tool calls -> typed
ChecklistItem list) with a monkeypatched `run_agent` -- a real run costs
real money on real literature/tool round-trips (verified manually: ~$0.3-0.6
on Haiku for a full 3-category checklist, see PROGRESS.md), so it isn't run
automatically on every test invocation. The extraction logic itself is real
code worth testing precisely, independent of the LLM's actual behavior.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from unittest.mock import patch

from multiome_agent.agent.checklist_generator import ChecklistItem, generate_checklist


@dataclass
class _FakeResult:
    answer: str = "done"
    tool_calls: list = field(default_factory=list)
    estimated_cost_usd: float = 0.01
    turn_count: int = 3
    hit_max_turns: bool = False


def _fake_call(category, gene, claim, pmid, journal, year, present, is_error=False):
    return {
        "name": "record_checklist_item",
        "input": {
            "category": category, "gene_or_motif": gene, "claim": claim,
            "pmid": pmid, "journal": journal, "year": year, "confirmed_present_in_data": present,
        },
        "is_error": is_error,
    }


def test_generate_checklist_extracts_only_record_checklist_item_calls():
    fake_result = _FakeResult(tool_calls=[
        {"name": "search_pubmed", "input": {"query": "x"}, "is_error": False},
        _fake_call("rna_marker", "CD14", "monocyte marker", "123", "Nature", "2020", True),
        {"name": "tf_motif_correlation", "input": {"gene": "SPI1"}, "is_error": False},
        _fake_call("motif", "SPI1", "myeloid motif", "456", "Cell", "2019", True),
    ])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, result = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert len(items) == 2
    assert all(isinstance(i, ChecklistItem) for i in items)
    assert items[0].category == "rna_marker" and items[0].gene_or_motif == "CD14"
    assert items[1].category == "motif" and items[1].pmid == "456"
    assert result is fake_result


def test_generate_checklist_excludes_error_tool_calls():
    fake_result = _FakeResult(tool_calls=[
        _fake_call("rna_marker", "BROKEN", "bad call", "1", "J", "2020", False, is_error=True),
        _fake_call("rna_marker", "CD14", "good call", "2", "Nature", "2020", True),
    ])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, _ = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert len(items) == 1
    assert items[0].gene_or_motif == "CD14"


def test_generate_checklist_handles_zero_items_gracefully():
    fake_result = _FakeResult(tool_calls=[{"name": "search_pubmed", "input": {}, "is_error": False}])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, _ = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert items == []


def test_generate_checklist_dedups_same_gene_within_category():
    # Real bug from an actual run: recorded two "tf_motif_tracking" items
    # both about FOXA1 (two different literature framings of the same
    # gene) while burning $6/60 turns searching for a genuinely different
    # 3rd candidate. Dedup is deterministic, not relied on the model's own
    # discipline -- keep only the first occurrence per (category, gene).
    fake_result = _FakeResult(tool_calls=[
        _fake_call("tf_motif_tracking", "FOXA1", "claim A", "1", "J1", "2020", True),
        _fake_call("tf_motif_tracking", "FOXA1", "claim B (different framing)", "2", "J2", "2023", True),
        _fake_call("tf_motif_tracking", "GATA3", "claim C", "3", "J3", "2021", True),
    ])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, _ = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert len(items) == 2
    assert [i.gene_or_motif for i in items] == ["FOXA1", "GATA3"]
    assert items[0].pmid == "1"  # first occurrence kept, not the later duplicate


def test_generate_checklist_allows_same_gene_across_different_categories():
    # Cross-category reuse is legitimate: "FOXA1 motif is elevated" (motif
    # category) and "FOXA1 expression tracks its own motif" (tf_motif_tracking
    # category) are two distinct, non-redundant claims about the same gene.
    fake_result = _FakeResult(tool_calls=[
        _fake_call("motif", "FOXA1", "motif claim", "1", "J1", "2020", True),
        _fake_call("tf_motif_tracking", "FOXA1", "tracking claim", "2", "J2", "2023", True),
    ])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, _ = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert len(items) == 2
    assert {i.category for i in items} == {"motif", "tf_motif_tracking"}


def test_generate_checklist_dedup_is_case_insensitive():
    fake_result = _FakeResult(tool_calls=[
        _fake_call("rna_marker", "Foxa1", "claim A", "1", "J1", "2020", True),
        _fake_call("rna_marker", "FOXA1", "claim B", "2", "J2", "2021", True),
    ])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, _ = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert len(items) == 1
