"""Tests for the report-generator orchestration (agent/report_generator.py).

`_render_section*` are pure functions, tested directly. `generate_report`'s
full orchestration is tested with every dependency mocked -- it makes real
LLM calls and real fault-injection reruns through several other modules,
so a real invocation is validated manually (see PROGRESS.md), not on every
test run.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from unittest.mock import MagicMock, patch

import pytest

from multiome_agent.agent.checklist_generator import ChecklistItem
from multiome_agent.agent.novelty import JudgedFinding
from multiome_agent.agent.report_generator import ReportCost, _extract_identity_blurb, _render_section1, _render_section2, _render_section3, _render_section4, _render_section5, generate_report


def test_extract_identity_blurb_real_pattern_with_double_asterisk_before_colon():
    # Real pattern captured from an actual run's reasoning trail: "---\n\n##
    # Summary\n\n**Dataset Identity**: <statement>\n\n**Known-Biology
    # Checklist**..." -- the raw regression this fixes: dumping the WHOLE
    # answer (including a checklist restatement redundant with the parsed
    # `items` list) into the report instead of just the identity statement.
    answer = (
        "---\n\n## Summary\n\n**Dataset Identity**: This is a **multi-cell-line "
        "pooled cancer dataset** containing epithelial cancer cells.\n\n"
        "**Known-Biology Checklist** (complete, all three categories):\n\n1. blah"
    )
    assert _extract_identity_blurb(answer) == (
        "This is a **multi-cell-line pooled cancer dataset** containing epithelial cancer cells."
    )


def test_extract_identity_blurb_real_pattern_with_colon_before_double_asterisk():
    answer = (
        "---\n\n## Summary\n\n**Dataset Identity:** This is a **PBMC multiome dataset** "
        "containing T cells.\n\n**Known-Biology Checklist**"
    )
    assert _extract_identity_blurb(answer) == "This is a **PBMC multiome dataset** containing T cells."


def test_extract_identity_blurb_no_heading_at_all_uses_leading_prose():
    assert _extract_identity_blurb("This is clearly PBMC data based on canonical markers.") == (
        "This is clearly PBMC data based on canonical markers."
    )


def test_extract_identity_blurb_no_identity_line_falls_through_to_remaining_prose():
    # Both a leading rule AND a leading heading are stripped, leaving
    # whatever real prose remains -- better than an empty-string fallback
    # when there's genuinely something left to show.
    text = _extract_identity_blurb("---\n\n## Summary\n\nNo identity line here at all, just a table.")
    assert text == "No identity line here at all, just a table."


def test_extract_identity_blurb_falls_back_when_truly_nothing_extractable():
    assert "reasoning trail" in _extract_identity_blurb("---")


def test_extract_identity_blurb_truncates_at_line_boundary_not_mid_table_row():
    # Real bug from an actual shareseq-multi-cell-lines run: the identity answer led
    # with "## Identity" (no "## Summary" wrapper, so the regex-match
    # branch never fires) followed by an 8-row markdown table: a plain
    # text[:900] cut landed mid-row, shipping a broken table (a row with
    # only 2 of 8 cells and a stray "...") in the committed report.
    rows = "\n".join(f"| Line{i} | ACH-{i:06d} | {700 - i} | Lineage{i} |" for i in range(20))
    answer = (
        "## Identity\n\n**Intro sentence about the pooled-panel identity.**"
        "\n\n| Cell line | DepMap ID | Cells | Lineage |\n|---|---|---|---|\n" + rows +
        "\n\nSome trailing prose.\n\n## Checklist (3 / 3 / 3)\n\nmore content"
    )
    # max_chars=150 lands inside the table's data rows (intro ~65 chars +
    # header/separator ~65 chars leaves room for one row before the cut) --
    # exactly the failure mode the old plain text[:max_chars] hit for real.
    blurb = _extract_identity_blurb(answer, max_chars=150)
    for line in blurb.rstrip(".").split("\n"):
        assert line.count("|") in (0, 5), f"cut mid-row: {line!r}"


def test_report_cost_aggregates_and_renders():
    cost = ReportCost()
    cost.add("a", 0.1)
    cost.add("b", 0.2)
    cost.add("a", 0.05)  # accumulates within the same step
    assert cost.total == pytest.approx(0.35)
    assert "a: $0.1500" in cost.render()
    assert "Total: $0.3500" in cost.render()


def test_render_section1_includes_loader_decision_and_summary():
    decision = MagicMock(strategy="combined_single_file", reasoning="both feature types in one file")
    text = _render_section1("tenx-cell-ranger", "500 cells, ...", decision)
    assert "combined_single_file" in text
    assert "500 cells" in text


def test_render_section2_handles_empty_and_populated():
    assert "no genes were checked" in _render_section2([])
    checks = [{"gene": "CD14", "rna_cluster": "1", "matched_atac_cluster": "1", "gene_activity_confirms": True}]
    text = _render_section2(checks)
    assert "CD14" in text and "1/1 checked markers cross-validate" in text


def test_render_section3_groups_by_category_and_flags_missing():
    items = [ChecklistItem("rna_marker", "CD14", "monocyte marker", "123", "Nature", "2020", True)]
    text = _render_section3(items, "This looks like PBMC based on markers.")
    assert "CD14" in text
    assert "no grounded, data-present candidate" in text  # motif/tf_motif_tracking empty


def _fake_model_results(detected, diag, false_alarm, cost, answer="Looks fine."):
    scen = [
        {"fault_type": "clean_control", "description": "clean", "answer": answer,
         "detected": false_alarm, "diagnosis_matches": None, "false_alarm": false_alarm},
        {"fault_type": "injected_doublets", "description": "doublets", "answer": "Doublets found.",
         "detected": detected, "diagnosis_matches": diag, "false_alarm": False},
    ]
    return {"scenario_results": scen, "total_cost_usd": cost}


def test_render_section4_renders_comparison_table_and_narrative_model_answers():
    results_by_model = {
        "claude-haiku-4-5": _fake_model_results(True, True, False, 0.02),
        "claude-opus-5": _fake_model_results(True, True, False, 0.10, answer="Looks fine (opus)."),
    }
    text = _render_section4(results_by_model, "claude-opus-5")
    assert "claude-haiku-4-5" in text and "claude-opus-5" in text
    assert "1/1" in text  # detected/diag for the 1 real fault scenario
    assert "Looks fine (opus)." in text
    assert "Looks fine." not in text  # only the narrative model's answers are shown in full


def test_render_section5_no_findings_and_negative_control():
    text = _render_section5([], {"n_findings_proposed": 0}, ["limitation A"])
    assert "No novel findings were proposed" in text
    assert "0 finding(s)" in text
    assert "limitation A" in text


def test_render_section5_missing_verdict_is_labeled():
    j = JudgedFinding(finding="X", evidence="ev", confidence="low", verdict=None, judge_cost_usd=0.01)
    text = _render_section5([j], {"n_findings_proposed": 0}, [])
    assert "NO VERDICT RECORDED" in text


def test_render_section5_negative_control_lists_both_judged_and_unjudged_findings():
    positive = JudgedFinding(
        finding="spurious cross-modal correlation", evidence="rho=0.14, p=0.002",
        confidence="medium", verdict={"verdict": "struck_down", "likely_artifact": True}, judge_cost_usd=1.19,
        finding_type="positive_relationship",
    )
    concern = {"finding": "cross-modal cluster matching collapses onto one partner", "evidence": "4 clusters -> same ATAC partner", "confidence": "medium", "finding_type": "no_signal_or_concern"}
    negative_control = {
        "n_findings_proposed": 2, "judged_positive_findings": [positive], "no_signal_or_concern_findings": [concern],
    }
    text = _render_section5([], negative_control, [])
    assert "spurious cross-modal correlation" in text and "struck_down" in text
    assert "cross-modal cluster matching collapses" in text
    assert "not adversarially judged" in text
    assert "no_signal_or_concern" in text


def test_render_section5_partial_verdict_missing_key_does_not_crash():
    # Real bug from an actual Opus run: record_judger_verdict's "required"
    # schema fields aren't server-enforced -- the model called the tool
    # without a "verdict" key at all, crashing the report AFTER every
    # expensive upstream step had already run and been paid for.
    partial = {"likely_artifact": True, "artifact_reasoning": "r1"}  # no "verdict", no "already_known"/"novelty_reasoning"
    j = JudgedFinding(finding="X", evidence="ev", confidence="high", verdict=partial, judge_cost_usd=0.02)
    text = _render_section5([j], {"n_findings_proposed": 0}, [])
    assert "NO VERDICT RECORDED" in text
    assert "(not provided)" in text
    assert "r1" in text


@dataclass
class _FakeAgentResult:
    answer: str = "PBMC identified."
    tool_calls: list = field(default_factory=list)
    estimated_cost_usd: float = 0.01


@dataclass
class _FakeScenarioResult:
    fault_type: str = "clean_control"
    severity: float = 0.0
    description: str = "clean"
    answer: str = "looks fine"
    detected: bool = False
    diagnosis_matches: object = None
    false_alarm: bool = False
    cost: float = 0.01
    turns: int = 2
    grounding: dict = field(default_factory=lambda: {"all_grounded": True})
    reasoning_trail_path: str = "x"


def test_generate_report_orchestration_writes_file_with_all_mocked(tmp_path):
    fake_mdata = MagicMock()
    fake_decision = MagicMock(strategy="combined_single_file", reasoning="one file, both feature types", cost_usd=0.001)
    checklist_items = [ChecklistItem("rna_marker", "CD14", "monocyte marker", "1", "Nature", "2020", True)]
    checklist_result = _FakeAgentResult(answer="This is PBMC.")
    findings = [{"finding": "F1", "evidence": "ev1", "confidence": "low"}]
    novelty_result = _FakeAgentResult()
    judged = [JudgedFinding(finding="F1", evidence="ev1", confidence="low", verdict={"verdict": "struck_down", "likely_artifact": True, "artifact_reasoning": "r", "already_known": False, "novelty_reasoning": "r2"}, judge_cost_usd=0.02)]
    fake_model_eval_result = {
        "model": "claude-haiku-4-5",
        "scenario_results": [dict(vars(_FakeScenarioResult()))],
        "total_cost_usd": 0.01,
    }

    with patch("multiome_agent.agent.report_generator.load_fixed_core_via_agent_decision", return_value=(fake_mdata, fake_decision)), \
         patch("multiome_agent.agent.report_generator._dataset_summary_text", return_value="500 cells summary"), \
         patch("multiome_agent.core.condition_detection.detect_condition_groups", return_value={"is_control_only": True, "condition_columns_found": {}}), \
         patch("multiome_agent.agent.report_generator.generate_checklist", return_value=(checklist_items, checklist_result)), \
         patch("multiome_agent.agent.report_generator.propose_novel_findings", return_value=(findings, novelty_result)), \
         patch("multiome_agent.agent.report_generator.judge_novel_findings", return_value=judged), \
         patch("multiome_agent.agent.report_generator.check_novelty_negative_control", return_value={"n_findings_proposed": 0, "findings": [], "cost_usd": 0.01, "answer": "none"}), \
         patch("multiome_agent.eval.eval_harness.build_scenarios", return_value=[{"fault_type": "clean_control"}]), \
         patch("multiome_agent.eval.eval_harness.run_model_eval", return_value=fake_model_eval_result), \
         patch("multiome_agent.eval.fault_injection.shuffle_rna_atac_pairing", return_value=(fake_mdata, MagicMock())), \
         patch("multiome_agent.agent.qc_summary.fixed_core_summary", return_value={}), \
         patch("multiome_agent.agent.qc_summary.format_qc_summary", return_value="shuffled summary"):
        result = generate_report("tenx-cell-ranger", model="claude-haiku-4-5", out_path=tmp_path / "tenx-cell-ranger-report.md")

    assert (tmp_path / "tenx-cell-ranger-report.md").exists()
    assert "PBMC" in result["report_text"] or "This is PBMC" in result["report_text"]
    assert result["cost"].total > 0
