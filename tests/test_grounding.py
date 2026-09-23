"""Tests for the grounding checker, using constructed fake run results (no
real API calls needed -- this only inspects text, it doesn't call Claude).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from multiome_agent.eval.grounding import check_grounding


@dataclass
class _FakeResult:
    answer: str
    tool_calls: list = field(default_factory=list)


def test_fully_grounded_answer_passes():
    result = _FakeResult(
        answer="SPI1 correlates with its motif at rho=0.595 (PMID 42749843 supports this).",
        tool_calls=[
            {"name": "tf_motif_correlation", "result_summary": '{"gene": "SPI1", "spearman_rho": 0.5948659415096681}'},
            {"name": "search_pubmed", "result_summary": '[{"pmid": "42749843", "title": "USP22..."}]'},
        ],
    )
    report = check_grounding(result)
    assert report["all_grounded"] is True
    assert report["numbers"][0]["grounded"] is True
    assert report["pmids"][0]["grounded"] is True


def test_number_not_in_any_tool_result_is_flagged():
    result = _FakeResult(
        answer="The correlation is rho=0.999, a very strong result.",
        tool_calls=[
            {"name": "tf_motif_correlation", "result_summary": '{"spearman_rho": 0.595}'},
        ],
    )
    report = check_grounding(result)
    assert report["all_grounded"] is False
    assert any(n["cited"] == "0.999" and not n["grounded"] for n in report["numbers"])


def test_fabricated_pmid_is_flagged():
    result = _FakeResult(
        answer="This is well established (PMID 99999999).",
        tool_calls=[
            {"name": "search_pubmed", "result_summary": '[{"pmid": "12345678", "title": "..."}]'},
        ],
    )
    report = check_grounding(result)
    assert report["all_grounded"] is False
    assert report["pmids"] == [{"pmid": "99999999", "grounded": False}]


def test_rounded_number_still_counts_as_grounded():
    # 0.5948659415096681 rounds to 0.59 at 2 decimals -- the cited value's
    # own precision, not an arbitrarily "close enough" tolerance.
    result = _FakeResult(
        answer="The correlation was approximately 0.59.",
        tool_calls=[{"name": "x", "result_summary": '{"spearman_rho": 0.5948659415096681}'}],
    )
    report = check_grounding(result)
    assert report["numbers"][0]["grounded"] is True


def test_wrong_value_at_matching_precision_is_not_falsely_grounded():
    # Same number of decimal places as the true value, but genuinely wrong
    # -- must not be waved through by the precision-aware rounding logic.
    result = _FakeResult(
        answer="The correlation was 0.612 exactly.",
        tool_calls=[{"name": "x", "result_summary": '{"spearman_rho": 0.5948659415096681}'}],
    )
    report = check_grounding(result)
    assert report["numbers"][0]["grounded"] is False


def test_unicode_scientific_notation_is_grounded_as_one_number():
    # LLM answers often render "3.5e-49" as "3.5 × 10⁻⁴⁹" -- must be treated
    # as one number (~3.5e-49), not decomposed into a bare "3.5" that then
    # fails to match against the tool's tiny p-value.
    result = _FakeResult(
        answer="This is highly significant (p = 3.5 × 10⁻⁴⁹).",
        tool_calls=[{"name": "x", "result_summary": '{"pvalue": 3.5279961399705646e-49}'}],
    )
    report = check_grounding(result)
    assert len(report["numbers"]) == 1
    assert report["numbers"][0]["grounded"] is True


def test_pmid_seen_only_via_fetch_pubmed_abstracts_is_grounded():
    # fetch_pubmed_abstracts' real return shape keys the dict BY pmid
    # ({"123": {"title": ..., "abstract": ...}}) -- no literal "pmid" field
    # at all, unlike search_pubmed's list-of-records shape. A PMID the agent
    # only ever saw this way must still count as grounded.
    result = _FakeResult(
        answer="This is supported by the literature (PMID 36662812).",
        tool_calls=[
            {
                "name": "fetch_pubmed_abstracts",
                "result": '{"36662812": {"title": "SWI/SNF Blockade...", "abstract": "In AML..."}}',
            },
        ],
    )
    report = check_grounding(result)
    assert report["all_grounded"] is True
    assert report["pmids"] == [{"pmid": "36662812", "grounded": True}]


def test_grounding_uses_full_result_field_not_truncated_summary():
    # A long tool result (e.g. several fetched abstracts) whose real PMID
    # would be truncated away in "result_summary" but is present in the
    # full "result" field -- grounding must read the latter.
    long_abstract = "x" * 2000
    full_result = f'{{"11111111": {{"abstract": "irrelevant"}}, "99999999": {{"abstract": "{long_abstract}"}}}}'
    result = _FakeResult(
        answer="See PMID 99999999 for details.",
        tool_calls=[
            {"name": "fetch_pubmed_abstracts", "result": full_result, "result_summary": full_result[:50] + "..."},
        ],
    )
    report = check_grounding(result)
    assert report["pmids"] == [{"pmid": "99999999", "grounded": True}]


def test_no_numbers_or_pmids_is_vacuously_grounded():
    result = _FakeResult(answer="This dataset shows expected PBMC biology.", tool_calls=[])
    report = check_grounding(result)
    assert report["all_grounded"] is True
    assert report["numbers"] == []
    assert report["pmids"] == []
