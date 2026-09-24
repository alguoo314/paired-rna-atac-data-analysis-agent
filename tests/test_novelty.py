"""Tests for novel-finding proposal + adversarial judging (agent/novelty.py).
Unit-tests the extraction logic with a monkeypatched `run_agent` -- real
runs cost real money on real tool/literature round-trips (see PROGRESS.md
for the manually-verified real end-to-end run).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from unittest.mock import patch

from multiome_agent.agent.novelty import check_novelty_negative_control, judge_novel_findings, propose_novel_findings


@dataclass
class _FakeResult:
    answer: str = "done"
    tool_calls: list = field(default_factory=list)
    estimated_cost_usd: float = 0.01
    turn_count: int = 3
    hit_max_turns: bool = False


def _finding_call(finding, evidence, confidence, finding_type="positive_relationship", is_error=False):
    return {
        "name": "record_novel_finding",
        "input": {"finding": finding, "evidence": evidence, "confidence": confidence, "finding_type": finding_type},
        "is_error": is_error,
    }


def _verdict_call(verdict, likely_artifact=False, already_known=False, is_error=False):
    return {
        "name": "record_judger_verdict",
        "input": {
            "likely_artifact": likely_artifact, "artifact_reasoning": "r1",
            "already_known": already_known, "novelty_reasoning": "r2", "verdict": verdict,
        },
        "is_error": is_error,
    }


def test_propose_novel_findings_extracts_only_record_calls():
    fake = _FakeResult(tool_calls=[
        {"name": "search_pubmed", "input": {}, "is_error": False},
        _finding_call("X correlates with Y", "rho=0.6 in cluster 3", "medium"),
    ])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        findings, result = propose_novel_findings(mdata=object(), qc_summary="qc", dataset_context="ctx")
    assert len(findings) == 1
    assert findings[0]["finding"] == "X correlates with Y"
    assert result is fake


def test_propose_novel_findings_zero_is_fine():
    fake = _FakeResult(tool_calls=[{"name": "search_pubmed", "input": {}, "is_error": False}])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        findings, _ = propose_novel_findings(mdata=object(), qc_summary="qc", dataset_context="ctx")
    assert findings == []


def test_propose_novel_findings_caps_positive_relationship_at_three():
    # Real run motivation: the model is asked to try up to 6 candidates but
    # stop recording after 3 that survive self-rejection -- don't trust that
    # discipline alone, cap deterministically. no_signal_or_concern findings
    # aren't judged downstream, so they're never capped.
    fake = _FakeResult(tool_calls=[
        _finding_call("F1", "ev1", "high"),
        _finding_call("F2", "ev2", "high"),
        _finding_call("F3", "ev3", "medium"),
        _finding_call("F4", "ev4", "low"),
        _finding_call("F5", "ev5", "low", finding_type="no_signal_or_concern"),
    ])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        findings, _ = propose_novel_findings(mdata=object(), qc_summary="qc", dataset_context="ctx")
    positive = [f for f in findings if f["finding_type"] == "positive_relationship"]
    assert len(positive) == 3
    assert [f["finding"] for f in positive] == ["F1", "F2", "F3"]
    assert any(f["finding"] == "F5" for f in findings)


def test_judge_novel_findings_pairs_verdict_with_original_finding():
    findings = [{"finding": "A", "evidence": "ev-A", "confidence": "high"}]
    fake = _FakeResult(tool_calls=[_verdict_call("struck_down", already_known=True)])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        judged = judge_novel_findings(findings, mdata=object(), qc_summary="qc", dataset_context="ctx")
    assert len(judged) == 1
    assert judged[0].finding == "A"
    assert judged[0].verdict["verdict"] == "struck_down"
    assert judged[0].verdict["already_known"] is True


def test_judge_novel_findings_handles_missing_verdict():
    findings = [{"finding": "A", "evidence": "ev-A", "confidence": "low"}]
    fake = _FakeResult(tool_calls=[])  # judger never called record_judger_verdict
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        judged = judge_novel_findings(findings, mdata=object(), qc_summary="qc", dataset_context="ctx")
    assert judged[0].verdict is None


def test_judge_novel_findings_backfills_missing_verdict_field_survives():
    # Real bug from an actual run: the judger called record_judger_verdict
    # with thorough, correct reasoning (likely_artifact=False,
    # already_known=False) but omitted the "verdict" enum field itself --
    # "required" isn't server-enforced. Backfill from the two booleans
    # rather than rendering "NO VERDICT RECORDED" for a real, reasoned call.
    findings = [{"finding": "A", "evidence": "ev-A", "confidence": "medium"}]
    fake = _FakeResult(tool_calls=[{
        "name": "record_judger_verdict",
        "input": {"likely_artifact": False, "artifact_reasoning": "r1", "already_known": False, "novelty_reasoning": "r2"},
        "is_error": False,
    }])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        judged = judge_novel_findings(findings, mdata=object(), qc_summary="qc", dataset_context="ctx")
    assert judged[0].verdict["verdict"] == "survives"


def test_judge_novel_findings_backfills_missing_verdict_field_struck_down():
    findings = [{"finding": "A", "evidence": "ev-A", "confidence": "medium"}]
    fake = _FakeResult(tool_calls=[{
        "name": "record_judger_verdict",
        "input": {"likely_artifact": True, "artifact_reasoning": "r1", "already_known": False, "novelty_reasoning": "r2"},
        "is_error": False,
    }])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        judged = judge_novel_findings(findings, mdata=object(), qc_summary="qc", dataset_context="ctx")
    assert judged[0].verdict["verdict"] == "struck_down"


def test_negative_control_reports_hallucination_when_findings_proposed():
    # No finding_type set -> defaults to positive_relationship -> gets
    # judged too; run_agent is mocked to return the same propose-shaped
    # fake for the judge call as well, so verdict ends up None (no
    # record_judger_verdict tool call in it) -- fine, this test only
    # checks the proposed count, not the judged verdict.
    fake = _FakeResult(tool_calls=[_finding_call("spurious link", "weak rho", "low")])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        result = check_novelty_negative_control(shuffled_mdata=object(), shuffled_qc_summary="qc", dataset_context="ctx")
    assert result["n_findings_proposed"] == 1


def test_negative_control_only_judges_positive_relationship_findings():
    # Real design: a negative-control run can propose a mix of a genuine
    # hallucination candidate (needs adversarial judging) and a correct
    # "nothing here" or tooling-concern observation (doesn't). Only the
    # former should trigger a Judger call.
    propose_fake = _FakeResult(tool_calls=[
        _finding_call("KLF2/KLF4 motif sign-split", "rho=+0.14/-0.11", "medium", finding_type="positive_relationship"),
        _finding_call("no coupling found for any lineage TF, as expected under shuffling", "all |rho|<0.09", "low", finding_type="no_signal_or_concern"),
    ])
    judge_fake = _FakeResult(tool_calls=[_verdict_call("struck_down", likely_artifact=True)])

    def fake_run_agent(question, **kwargs):
        # "skeptical" alone isn't a safe discriminator -- the propose
        # question's own negative-control guidance now mentions "a
        # skeptical judge" in passing. Match the judge question's unique
        # instruction to call the verdict tool instead.
        return judge_fake if "record_judger_verdict` exactly once" in question else propose_fake

    with patch("multiome_agent.agent.novelty.run_agent", side_effect=fake_run_agent):
        result = check_novelty_negative_control(shuffled_mdata=object(), shuffled_qc_summary="qc", dataset_context="ctx")

    assert result["n_findings_proposed"] == 2
    assert len(result["judged_positive_findings"]) == 1
    assert result["judged_positive_findings"][0].verdict["verdict"] == "struck_down"
    assert len(result["no_signal_or_concern_findings"]) == 1
    assert result["no_signal_or_concern_findings"][0]["finding"] == "no coupling found for any lineage TF, as expected under shuffling"


def test_negative_control_skips_judging_entirely_when_all_concerns():
    propose_fake = _FakeResult(tool_calls=[
        _finding_call("cross-modal cluster matching is degenerate", "4 clusters -> 1 ATAC partner", "medium", finding_type="no_signal_or_concern"),
    ])
    calls = []

    def fake_run_agent(question, **kwargs):
        calls.append(question)
        return propose_fake

    with patch("multiome_agent.agent.novelty.run_agent", side_effect=fake_run_agent):
        result = check_novelty_negative_control(shuffled_mdata=object(), shuffled_qc_summary="qc", dataset_context="ctx")

    assert len(calls) == 1  # only the propose call -- no Judger call spent on a non-positive finding
    assert result["judged_positive_findings"] == []
    assert len(result["no_signal_or_concern_findings"]) == 1


def test_propose_novel_findings_defaults_missing_finding_type_to_positive_relationship():
    # "required" isn't server-enforced -- a model could omit finding_type.
    fake = _FakeResult(tool_calls=[
        {"name": "record_novel_finding", "input": {"finding": "X", "evidence": "Y", "confidence": "low"}, "is_error": False},
    ])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        findings, _ = propose_novel_findings(mdata=object(), qc_summary="qc", dataset_context="ctx")
    assert findings[0]["finding_type"] == "positive_relationship"


def test_negative_control_clean_when_no_findings_proposed():
    fake = _FakeResult(tool_calls=[])
    with patch("multiome_agent.agent.novelty.run_agent", return_value=fake):
        result = check_novelty_negative_control(shuffled_mdata=object(), shuffled_qc_summary="qc", dataset_context="ctx")
    assert result["n_findings_proposed"] == 0
