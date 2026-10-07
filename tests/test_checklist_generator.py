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

from multiome_agent.agent.checklist_generator import CHECKLIST_CATEGORIES, ChecklistItem, generate_checklist


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
    # Real bug from an actual run: recorded two "peak_to_gene" items both
    # about FOXA1 (two different literature framings of the same gene)
    # while burning $6/60 turns searching for a genuinely different 3rd
    # candidate. Dedup is deterministic, not relied on the model's own
    # discipline -- keep only the first occurrence per (category, gene).
    fake_result = _FakeResult(tool_calls=[
        _fake_call("peak_to_gene", "FOXA1", "claim A", "1", "J1", "2020", True),
        _fake_call("peak_to_gene", "FOXA1", "claim B (different framing)", "2", "J2", "2023", True),
        _fake_call("peak_to_gene", "GATA3", "claim C", "3", "J3", "2021", True),
    ])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, _ = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert len(items) == 2
    assert [i.gene_or_motif for i in items] == ["FOXA1", "GATA3"]
    assert items[0].pmid == "1"  # first occurrence kept, not the later duplicate


def test_generate_checklist_allows_same_gene_across_different_categories():
    # Cross-category reuse is legitimate: "FOXA1 motif is elevated" (motif
    # category) and "FOXA1 regulates a distal target" (peak_to_gene
    # category) are two distinct, non-redundant claims about the same gene.
    fake_result = _FakeResult(tool_calls=[
        _fake_call("motif", "FOXA1", "motif claim", "1", "J1", "2020", True),
        _fake_call("peak_to_gene", "FOXA1", "distal-link claim", "2", "J2", "2023", True),
    ])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, _ = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert len(items) == 2
    assert {i.category for i in items} == {"motif", "peak_to_gene"}


def test_generate_checklist_dedup_is_case_insensitive():
    fake_result = _FakeResult(tool_calls=[
        _fake_call("rna_marker", "Foxa1", "claim A", "1", "J1", "2020", True),
        _fake_call("rna_marker", "FOXA1", "claim B", "2", "J2", "2021", True),
    ])
    with patch("multiome_agent.agent.checklist_generator.run_agent", return_value=fake_result):
        items, _ = generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    assert len(items) == 1


def test_generate_checklist_shortlist_fully_covers_category_skips_agent_entirely():
    # When the cheap shortlist path reaches 3 confirmed items in a category,
    # the expensive run_agent must not be called for it at all.
    shortlist_item = ChecklistItem("rna_marker", "CD14", "claim", "1", "J", "2020", True)
    with patch("multiome_agent.agent.checklist_generator.resolve_cheap_identity", return_value=["T-47D"]), \
         patch("multiome_agent.agent.checklist_generator.shortlist_checklist_items_for_category",
               return_value=[shortlist_item, shortlist_item, shortlist_item]) as mock_shortlist, \
         patch("multiome_agent.agent.checklist_generator.run_agent") as mock_run_agent:
        items, result = generate_checklist(
            mdata=object(), qc_summary="qc", dataset_context="ctx", categories=("rna_marker",),
        )

    mock_run_agent.assert_not_called()
    assert mock_shortlist.call_count == 1
    assert len(items) == 1  # 3 identical items dedup to 1 (same category, same gene)
    assert result.estimated_cost_usd == 0.0


def test_generate_checklist_shortlist_shortfall_triggers_scoped_fallback():
    # Shortlist finds only 1/3 for "motif" -- run_agent must be called, and
    # ONLY for the category that still needs filling.
    with patch("multiome_agent.agent.checklist_generator.resolve_cheap_identity", return_value=["T-47D"]), \
         patch("multiome_agent.agent.checklist_generator.shortlist_checklist_items_for_category",
               return_value=[ChecklistItem("motif", "FOXA1", "claim", "1", "J", "2020", True)]), \
         patch("multiome_agent.agent.checklist_generator.run_agent") as mock_run_agent:
        mock_run_agent.return_value = _FakeResult(tool_calls=[
            _fake_call("motif", "GATA3", "claim B", "2", "J2", "2021", True),
        ])
        items, result = generate_checklist(
            mdata=object(), qc_summary="qc", dataset_context="ctx", categories=("motif",),
        )

    mock_run_agent.assert_called_once()
    fallback_question = mock_run_agent.call_args[0][0]
    assert "motif" in fallback_question
    assert "rna_marker" not in fallback_question  # only the categories still needing work are asked about
    assert "T-47D" in fallback_question  # resolved identity passed through, not re-derived
    assert {i.gene_or_motif for i in items} == {"FOXA1", "GATA3"}


def test_generate_checklist_no_identity_available_falls_back_to_full_agent_unconditionally():
    with patch("multiome_agent.agent.checklist_generator.resolve_cheap_identity", return_value=None), \
         patch("multiome_agent.agent.checklist_generator.shortlist_checklist_items_for_category") as mock_shortlist, \
         patch("multiome_agent.agent.checklist_generator.run_agent") as mock_run_agent:
        mock_run_agent.return_value = _FakeResult(tool_calls=[])
        generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx")

    mock_shortlist.assert_not_called()
    mock_run_agent.assert_called_once()


def test_generate_checklist_use_shortlist_false_skips_shortlist_even_with_identity():
    with patch("multiome_agent.agent.checklist_generator.resolve_cheap_identity") as mock_resolve, \
         patch("multiome_agent.agent.checklist_generator.run_agent") as mock_run_agent:
        mock_run_agent.return_value = _FakeResult(tool_calls=[])
        generate_checklist(mdata=object(), qc_summary="qc", dataset_context="ctx", use_shortlist=False)

    mock_resolve.assert_not_called()
    mock_run_agent.assert_called_once()


def test_generate_checklist_existing_items_preserved_for_untouched_categories():
    # Categories outside `categories` must pass through completely untouched
    # -- no shortlist call, no agent call, same objects back out.
    old_marker_item = ChecklistItem("rna_marker", "CD14", "old claim", "1", "J", "2020", True)
    old_motif_item = ChecklistItem("motif", "SPI1", "old claim", "2", "J", "2020", True)
    with patch("multiome_agent.agent.checklist_generator.resolve_cheap_identity", return_value=["T-47D"]), \
         patch("multiome_agent.agent.checklist_generator.shortlist_checklist_items_for_category",
               return_value=[ChecklistItem("regulon_target", "SPI1->CD14", "new claim", "3", "J", "2021", True)] * 3), \
         patch("multiome_agent.agent.checklist_generator.run_agent") as mock_run_agent:
        items, _ = generate_checklist(
            mdata=object(), qc_summary="qc", dataset_context="ctx",
            categories=("regulon_target",), existing_items=[old_marker_item, old_motif_item],
        )

    mock_run_agent.assert_not_called()
    assert old_marker_item in items
    assert old_motif_item in items
    assert any(i.category == "regulon_target" for i in items)
    assert len(items) == 3  # 2 preserved + 1 deduped-from-3 new


def test_generate_checklist_default_categories_cover_all_four():
    # "tf_motif_tracking" (TF-vs-own-motif) was deliberately removed from the
    # checklist in favor of "regulon_target" (TF-vs-specific-other-gene) --
    # a structurally different, harder-to-already-know claim. The tool
    # (tf_motif_correlation) remains real and usable elsewhere; it's just
    # no longer a checklist category.
    assert CHECKLIST_CATEGORIES == ("rna_marker", "motif", "peak_to_gene", "regulon_target")
