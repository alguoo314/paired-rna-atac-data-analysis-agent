"""Tests for the drug-comparison agent orchestration (agent/drug_comparison.py).

Real `run_agent`/`propose_novel_findings` calls cost real money, so they're
mocked here -- same convention as test_checklist_generator.py/test_novelty.py.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from unittest.mock import patch

from multiome_agent.agent.drug_comparison import _drug_comparison_question, propose_drug_comparison_findings


@dataclass
class _FakeResult:
    answer: str = "done"
    tool_calls: list = field(default_factory=list)
    estimated_cost_usd: float = 0.01


def test_propose_drug_comparison_findings_is_a_zero_cost_noop_when_control_only():
    condition_check = {"is_control_only": True, "condition_columns_found": {}}
    with patch("multiome_agent.agent.drug_comparison.propose_novel_findings") as mock_propose:
        findings, result = propose_drug_comparison_findings(mdata=object(), condition_check=condition_check, qc_summary="qc")

    mock_propose.assert_not_called()  # no agent call, no API cost, at all
    assert findings == []
    assert result.estimated_cost_usd == 0.0
    assert result.turn_count == 0


def test_propose_drug_comparison_findings_delegates_with_custom_question_when_real_axis_found():
    condition_check = {
        "is_control_only": False,
        "condition_columns_found": {"rna.drug_treatment": {"DMSO": 60, "Drug_X": 40}},
    }
    fake_findings = [{"finding": "F1", "evidence": "ev1", "confidence": "medium", "finding_type": "positive_relationship"}]
    fake_result = _FakeResult()
    with patch("multiome_agent.agent.drug_comparison.propose_novel_findings", return_value=(fake_findings, fake_result)) as mock_propose:
        findings, result = propose_drug_comparison_findings(mdata=object(), condition_check=condition_check, qc_summary="qc")

    mock_propose.assert_called_once()
    _, kwargs = mock_propose.call_args
    assert "custom_question" in kwargs
    assert "rna.drug_treatment" in kwargs["custom_question"]
    assert "check_batch_confound" in kwargs["custom_question"]
    assert findings is fake_findings
    assert result is fake_result


def test_drug_comparison_question_names_the_real_columns_and_group_sizes():
    condition_check = {"condition_columns_found": {"rna.drug_treatment": {"DMSO": 60, "Drug_X": 40}}}
    question = _drug_comparison_question(condition_check)
    assert "rna.drug_treatment" in question
    assert "DMSO" in question and "Drug_X" in question
    assert "cell_line" in question  # per-cell-line scoping instruction present
    assert "check_batch_confound" in question
