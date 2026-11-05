"""Tests for the QC-verdict step (agent/qc_verdict.py).

Unit-tests the extraction logic (record_qc_verdict tool call -> typed
QCVerdict) with a monkeypatched `run_agent` -- a real run costs real money
on real tool round-trips, so it isn't run automatically on every test
invocation (same convention as test_checklist_generator.py/test_novelty.py).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from unittest.mock import patch

from multiome_agent.agent.qc_verdict import QC_VERDICT_TOOLS, QCVerdict, generate_qc_verdict


@dataclass
class _FakeResult:
    answer: str = "done"
    tool_calls: list = field(default_factory=list)
    estimated_cost_usd: float = 0.01


def _verdict_call(status, reasons, is_error=False):
    return {"name": "record_qc_verdict", "input": {"status": status, "reasons": reasons}, "is_error": is_error}


def test_generate_qc_verdict_extracts_status_and_reasons():
    fake_result = _FakeResult(tool_calls=[
        {"name": "get_qc_summary", "input": {}, "is_error": False},
        _verdict_call("usable with caveats", ["elevated doublet rate", "thin A-673 arm"]),
    ])
    with patch("multiome_agent.agent.qc_verdict.run_agent", return_value=fake_result):
        verdict, result = generate_qc_verdict(mdata=object(), qc_summary="qc")

    assert isinstance(verdict, QCVerdict)
    assert verdict.status == "usable with caveats"
    assert verdict.reasons == ["elevated doublet rate", "thin A-673 arm"]
    assert result is fake_result


def test_generate_qc_verdict_excludes_error_tool_calls_takes_last_valid():
    fake_result = _FakeResult(tool_calls=[
        _verdict_call("do not use", ["bad call"], is_error=True),
        _verdict_call("usable", ["clean numbers"]),
    ])
    with patch("multiome_agent.agent.qc_verdict.run_agent", return_value=fake_result):
        verdict, _ = generate_qc_verdict(mdata=object(), qc_summary="qc")
    assert verdict.status == "usable"


def test_generate_qc_verdict_defaults_to_cautious_status_when_never_called():
    fake_result = _FakeResult(tool_calls=[{"name": "get_qc_summary", "input": {}, "is_error": False}])
    with patch("multiome_agent.agent.qc_verdict.run_agent", return_value=fake_result):
        verdict, _ = generate_qc_verdict(mdata=object(), qc_summary="qc")
    assert verdict.status == "usable with caveats"
    assert len(verdict.reasons) >= 1


def test_qc_verdict_tools_is_scoped_to_marker_and_cross_modal_tools_only():
    names = {t["name"] for t in QC_VERDICT_TOOLS}
    assert names == {
        "get_qc_summary", "list_clusters", "top_cluster_markers",
        "top_gene_activity_markers", "cross_modal_marker_check", "record_qc_verdict",
    }
    # Deliberately excludes the heavier analysis-menu/literature tools -- this
    # gate answers "is the data usable," not "what novel biology is in it."
    assert "regulon_inference" not in names
    assert "search_pubmed" not in names


def test_generate_qc_verdict_passes_require_tool_call(monkeypatch):
    captured = {}

    def _fake_run_agent(question, **kwargs):
        captured.update(kwargs)
        return _FakeResult(tool_calls=[_verdict_call("usable", ["fine"])])

    monkeypatch.setattr("multiome_agent.agent.qc_verdict.run_agent", _fake_run_agent)
    generate_qc_verdict(mdata=object(), qc_summary="qc")
    assert captured["require_tool_call"] == "record_qc_verdict"
    assert captured["tools"] is QC_VERDICT_TOOLS
