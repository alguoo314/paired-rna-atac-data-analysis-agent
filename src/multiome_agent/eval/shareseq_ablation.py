"""Ablation study for the shareseq-multi-cell-lines data: fixed-core-only vs.
core + agent-chosen analyses. CLAUDE.md: "Does agent choice actually help?"
-- answered empirically here, not assumed either way.
"""

from __future__ import annotations

from multiome_agent.agent.loop import run_agent
from multiome_agent.agent.shareseq_qc_summary import format_shareseq_qc_summary, shareseq_fixed_core_summary
from multiome_agent.agent.prompts import OWN_DATA_CONTEXT

QUESTION = (
    "This is a pooled multi-cell-line single-cell multiome dataset. Using the QC summary "
    "tool, assess data quality and clustering structure. Does the clustering look like it's "
    "capturing real, distinct biological populations, or could this be a technical artifact? "
    "If there's a specific angle worth investigating further, use the available tools."
)


def fixed_core_only_report(mdata) -> str:
    """What a reader sees with NO agent involvement -- the deterministic
    numbers alone, no interpretation, no literature context, no choice of
    what to look into next."""
    summary = shareseq_fixed_core_summary(mdata)
    return format_shareseq_qc_summary(summary)


def run_ablation(mdata, model: str | None = None) -> dict:
    """Runs both conditions for real and returns both, for direct comparison."""
    fixed_core_text = fixed_core_only_report(mdata)
    summary = shareseq_fixed_core_summary(mdata)
    qc_summary_text = format_shareseq_qc_summary(summary)

    agent_result = run_agent(
        QUESTION, model=model, mdata=mdata, qc_summary=qc_summary_text, dataset_context=OWN_DATA_CONTEXT
    )

    return {
        "fixed_core_only": fixed_core_text,
        "agent_answer": agent_result.answer,
        "agent_tool_calls": [tc["name"] for tc in agent_result.tool_calls],
        "agent_cost_usd": agent_result.estimated_cost_usd,
        "agent_turns": agent_result.turn_count,
    }
