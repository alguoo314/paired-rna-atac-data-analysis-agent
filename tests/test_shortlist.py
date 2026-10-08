"""Tests for the cost-saving shortlist-first path (agent/shortlist.py).
Every LLM call (`_call_shortlist_model`) and every network call
(`search_pubmed`/`fetch_pubmed_abstracts`) is mocked -- none of this spends
real API cost, matching this project's existing convention for testing the
expensive `run_agent` path (see test_checklist_generator.py/test_novelty.py).
"""

from __future__ import annotations

from unittest.mock import patch

from multiome_agent.agent.shortlist import (
    cheap_already_known_check,
    generate_shortlist_candidates,
    judge_abstract_supports_claim,
    resolve_cheap_identity,
    shortlist_checklist_items_for_category,
    verify_candidate_in_data,
)


class _FakeMuData:
    def __init__(self, found: dict):
        self._found = found
        self.mod = {}

    def mod_get(self):
        return self.mod


def test_resolve_cheap_identity_prefers_human_readable_cell_line_name(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.detect_identity_columns",
        lambda mdata: {"identity_columns_found": {
            "rna.cell_lines": {"ACH-000416": 1626, "ACH-000667": 817},
            "rna.cell_line_name": {"NCI-H838": 1626, "HCC-44": 817},
        }},
    )
    result = resolve_cheap_identity(object())
    assert result == ["NCI-H838", "HCC-44"]  # most abundant first, ACH-coded column ignored


def test_resolve_cheap_identity_falls_back_to_generic_non_numeric_column(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.detect_identity_columns",
        lambda mdata: {"identity_columns_found": {
            "rna.may_have_wrong_cell_line_label_based_on_rna_cluster": {"no": 100},
            "atac.donor_id": {"DonorA": 50, "DonorB": 30},
        }},
    )
    result = resolve_cheap_identity(object())
    assert result == ["DonorA", "DonorB"]


def test_resolve_cheap_identity_returns_none_when_nothing_found(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.detect_identity_columns",
        lambda mdata: {"identity_columns_found": {}},
    )
    assert resolve_cheap_identity(object()) is None


def test_resolve_cheap_identity_is_defensive_on_bad_mdata():
    # No monkeypatch -- the real detect_identity_columns will raise on a bare
    # object() (no `.mod`). Must return None, never propagate the exception --
    # this is a cost optimization, it must never be why the pipeline crashes.
    assert resolve_cheap_identity(object()) is None


def test_generate_shortlist_candidates_truncates_to_n(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist._call_shortlist_model",
        lambda prompt, tool, model: {"candidates": [
            {"primary": "A", "claim": "claim A"}, {"primary": "B", "claim": "claim B"},
            {"primary": "C", "claim": "claim C"},
        ]},
    )
    result = generate_shortlist_candidates("T-47D", "rna_marker", n=2)
    assert len(result) == 2
    assert result[0]["primary"] == "A"


def test_generate_shortlist_candidates_handles_none_result(monkeypatch):
    monkeypatch.setattr("multiome_agent.agent.shortlist._call_shortlist_model", lambda prompt, tool, model: None)
    assert generate_shortlist_candidates("T-47D", "rna_marker") == []


def test_judge_abstract_supports_claim_true_and_false(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist._call_shortlist_model",
        lambda prompt, tool, model: {"supports_claim": True},
    )
    assert judge_abstract_supports_claim("claim", "abstract") is True

    monkeypatch.setattr(
        "multiome_agent.agent.shortlist._call_shortlist_model",
        lambda prompt, tool, model: {"supports_claim": False},
    )
    assert judge_abstract_supports_claim("claim", "abstract") is False


def test_verify_candidate_in_data_dispatches_rna_marker(monkeypatch):
    monkeypatch.setattr("multiome_agent.agent.shortlist.marker_is_recovered", lambda rna, gene: gene == "CD14")
    mdata = type("M", (), {"mod": {"rna": object()}})()
    assert verify_candidate_in_data(mdata, "rna_marker", "CD14") is True
    assert verify_candidate_in_data(mdata, "rna_marker", "XYZ") is False


def test_verify_candidate_in_data_dispatches_regulon_target(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.infer_regulon_targets",
        lambda mdata, tf, target_gene=None, cell_line=None: {"targets": [{"gene": "CD14", "significant": True}, {"gene": "LYN", "significant": False}]},
    )
    assert verify_candidate_in_data(object(), "regulon_target", "SPI1->CD14") is True
    assert verify_candidate_in_data(object(), "regulon_target", "SPI1->LYN") is False
    assert verify_candidate_in_data(object(), "regulon_target", "SPI1->NOTAREALGENE") is False


def test_verify_candidate_in_data_regulon_target_malformed_primary_is_false():
    assert verify_candidate_in_data(object(), "regulon_target", "SPI1_no_arrow") is False


def test_verify_candidate_in_data_exception_is_caught(monkeypatch):
    def _raise(*a, **k):
        raise RuntimeError("boom")
    monkeypatch.setattr("multiome_agent.agent.shortlist.marker_is_recovered", _raise)
    assert verify_candidate_in_data(object(), "rna_marker", "CD14") is False


def test_shortlist_checklist_items_for_category_full_pipeline(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.generate_shortlist_candidates",
        lambda identity, category, n, model: [
            {"primary": "GOOD", "claim": "good claim"},
            {"primary": "FAILS_DATA", "claim": "fails data claim"},
            {"primary": "FAILS_LIT", "claim": "fails literature claim"},
            {"primary": "FAILS_JUDGE", "claim": "fails judge claim"},
        ],
    )
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.verify_candidate_in_data",
        lambda mdata, category, primary, cell_line=None: primary != "FAILS_DATA",
    )
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.search_pubmed",
        lambda claim: [] if claim == "fails literature claim" else [{"pmid": "123", "title": "t"}],
    )
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.fetch_pubmed_abstracts",
        lambda pmids: {"123": {"abstract": "real abstract text", "journal": "J", "year": "2020"}},
    )
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.judge_abstract_supports_claim",
        lambda claim, abstract, model: claim != "fails judge claim",
    )

    items = shortlist_checklist_items_for_category(object(), ["T-47D"], "rna_marker", needed=3)

    assert len(items) == 1
    assert items[0].gene_or_motif == "GOOD"
    assert items[0].pmid == "123"
    assert items[0].confirmed_present_in_data is True


def test_shortlist_checklist_items_for_category_stops_once_needed_reached(monkeypatch):
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.generate_shortlist_candidates",
        lambda identity, category, n, model: [{"primary": f"G{i}", "claim": f"claim {i}"} for i in range(5)],
    )
    monkeypatch.setattr("multiome_agent.agent.shortlist.verify_candidate_in_data", lambda *a, **k: True)
    monkeypatch.setattr("multiome_agent.agent.shortlist.search_pubmed", lambda claim: [{"pmid": "1", "title": "t"}])
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.fetch_pubmed_abstracts",
        lambda pmids: {"1": {"abstract": "x", "journal": "J", "year": "2020"}},
    )
    monkeypatch.setattr("multiome_agent.agent.shortlist.judge_abstract_supports_claim", lambda *a, **k: True)

    items = shortlist_checklist_items_for_category(object(), ["T-47D"], "rna_marker", needed=2)
    assert len(items) == 2


def test_shortlist_checklist_items_for_category_rotates_lines_and_records_free_replication(monkeypatch):
    # Each line's candidate-generation call is tracked so we can assert LINE_C's
    # is never spent once `needed` is already reached by LINE_A + LINE_B's claims
    # -- the rotation stops early, it doesn't exhaustively poll every line.
    generation_calls = []

    def _fake_generate(identity, category, n, model):
        generation_calls.append(identity)
        return {
            "LINE_A": [{"primary": "GENE1", "claim": "claim1"}],
            "LINE_B": [{"primary": "GENE2", "claim": "claim2"}],
            "LINE_C": [{"primary": "GENE3", "claim": "claim3"}],
        }[identity]

    monkeypatch.setattr("multiome_agent.agent.shortlist.generate_shortlist_candidates", _fake_generate)

    def _fake_verify(mdata, category, primary, cell_line=None):
        # GENE1 (framed around LINE_A) also replicates for free in LINE_C but not
        # LINE_B; GENE2 (framed around LINE_B) replicates nowhere else.
        if primary == "GENE1":
            return cell_line in ("LINE_A", "LINE_C")
        if primary == "GENE2":
            return cell_line == "LINE_B"
        return False

    monkeypatch.setattr("multiome_agent.agent.shortlist.verify_candidate_in_data", _fake_verify)
    monkeypatch.setattr("multiome_agent.agent.shortlist.search_pubmed", lambda claim: [{"pmid": "1", "title": "t"}])
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.fetch_pubmed_abstracts",
        lambda pmids: {"1": {"abstract": "x", "journal": "J", "year": "2020"}},
    )
    monkeypatch.setattr("multiome_agent.agent.shortlist.judge_abstract_supports_claim", lambda *a, **k: True)

    items = shortlist_checklist_items_for_category(
        object(), ["LINE_A", "LINE_B", "LINE_C"], "peak_to_gene", needed=2,
    )

    assert [i.gene_or_motif for i in items] == ["GENE1", "GENE2"]
    assert items[0].cell_lines == ["LINE_A", "LINE_C"]  # origin line first, then free-checked replication
    assert items[1].cell_lines == ["LINE_B"]
    assert "LINE_C" not in generation_calls  # needed=2 reached before LINE_C's own candidates were ever generated


def test_cheap_already_known_check_resolves_when_abstract_supports(monkeypatch):
    monkeypatch.setattr("multiome_agent.agent.shortlist.search_pubmed", lambda q: [{"pmid": "99", "title": "t"}])
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.fetch_pubmed_abstracts",
        lambda pmids: {"99": {"abstract": "supports it", "journal": "J", "year": "2021"}},
    )
    monkeypatch.setattr("multiome_agent.agent.shortlist.judge_abstract_supports_claim", lambda *a, **k: True)

    verdict = cheap_already_known_check("some finding")
    assert verdict["verdict"] == "struck_down"
    assert verdict["already_known"] is True


def test_cheap_already_known_check_returns_none_when_no_hits(monkeypatch):
    monkeypatch.setattr("multiome_agent.agent.shortlist.search_pubmed", lambda q: [])
    assert cheap_already_known_check("some finding") is None


def test_cheap_already_known_check_returns_none_when_judge_disagrees(monkeypatch):
    monkeypatch.setattr("multiome_agent.agent.shortlist.search_pubmed", lambda q: [{"pmid": "99", "title": "t"}])
    monkeypatch.setattr(
        "multiome_agent.agent.shortlist.fetch_pubmed_abstracts",
        lambda pmids: {"99": {"abstract": "unrelated", "journal": "J", "year": "2021"}},
    )
    monkeypatch.setattr("multiome_agent.agent.shortlist.judge_abstract_supports_claim", lambda *a, **k: False)
    assert cheap_already_known_check("some finding") is None
