"""Tests for the shareseq-multi-cell-lines-data ablation study. Skips gracefully if the
shareseq-multi-cell-lines data isn't configured locally. The real end-to-end test makes one
real, cheap Haiku call (consistent with the rest of this project's real-API
test pattern) -- not mocked, since the point is verifying real agent
behavior against real data, not just that the code runs.
"""

from __future__ import annotations

import pytest

from multiome_agent.agent.shareseq_fixed_core_cache import get_shareseq_fixed_core
from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD
from multiome_agent.eval.shareseq_ablation import fixed_core_only_report, run_ablation

SHARESEQ_DATA_CONFIGURED = bool(SHARESEQ_RNA_H5AD and SHARESEQ_ATAC_H5AD and SHARESEQ_RNA_HVG_H5AD)

pytestmark = pytest.mark.skipif(
    not SHARESEQ_DATA_CONFIGURED,
    reason="config/local_paths.yaml not set up locally; shareseq-multi-cell-lines data unavailable",
)


@pytest.fixture(scope="module")
def shareseq_fixed_core_mdata():
    # Cached (see agent/shareseq_fixed_core_cache.py) -- the pipeline now
    # includes chromVAR motif deviations (~39 min), too slow to recompute
    # per test session.
    return get_shareseq_fixed_core()


def test_fixed_core_only_report_has_no_agent_involvement(shareseq_fixed_core_mdata):
    report = fixed_core_only_report(shareseq_fixed_core_mdata)
    assert "cell lines" in report
    assert "ARI" in report


def test_real_ablation_run_produces_both_conditions(shareseq_fixed_core_mdata):
    result = run_ablation(shareseq_fixed_core_mdata, model="claude-haiku-4-5")
    assert result["fixed_core_only"]
    assert result["agent_answer"].strip() != ""
    assert result["agent_cost_usd"] > 0
    assert "get_qc_summary" in result["agent_tool_calls"]
