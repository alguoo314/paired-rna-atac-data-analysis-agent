"""QC-verdict step: after the fixed core is computed, run ONE agent call
whose only job is to decide whether this run's QC/clustering results are
usable at all, before any of the more expensive known-biology/novelty steps
spend money assuming they are.

Deliberately scoped to a SMALL tool subset (the QC summary itself plus the
marker/cross-modal-validation tools) rather than the full analysis-menu --
answering "is this data usable" needs cluster recoverability and marker
sanity-checks, not literature RAG or regulon/peak-to-gene correlations,
which are downstream steps that presuppose the data already passed this
gate. Keeping this agent small keeps the gate cheap.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from multiome_agent.agent.loop import TOOLS, AgentRunResult, run_agent
from multiome_agent.agent.prompts import OWN_DATA_CONTEXT
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

QC_VERDICT_STATUSES = ("usable", "usable with caveats", "do not use")

_QC_VERDICT_TOOL_NAMES = {
    "get_qc_summary", "list_clusters", "top_cluster_markers",
    "top_gene_activity_markers", "cross_modal_marker_check",
}

RECORD_QC_VERDICT_TOOL = {
    "name": "record_qc_verdict",
    "description": (
        "Record your final verdict on whether this dataset's fixed-core QC and clustering "
        "results are usable for downstream analysis. Call this exactly once, after checking "
        "real numbers (get_qc_summary) and real cluster structure (list_clusters, "
        "top_cluster_markers/top_gene_activity_markers, cross_modal_marker_check) -- not a "
        "guess from the QC numbers alone."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "status": {
                "type": "string",
                "enum": list(QC_VERDICT_STATUSES),
                "description": (
                    "'usable': QC numbers and cluster structure both look sound, no real "
                    "concerns. 'usable with caveats': workable overall, but with specific "
                    "named caveats worth a reader's attention (e.g. an elevated doublet rate, "
                    "a thin cell-line/cluster arm, a marker that doesn't cross-validate). "
                    "'do not use': a fundamental problem -- e.g. clustering fails to recover "
                    "any real structure, QC numbers are physically implausible, or cross-modal "
                    "agreement has collapsed."
                ),
            },
            "reasons": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Specific, numbers-grounded reasons for this verdict -- cite real tool results, don't guess.",
            },
        },
        "required": ["status", "reasons"],
    },
}

QC_VERDICT_TOOLS = [t for t in TOOLS if t["name"] in _QC_VERDICT_TOOL_NAMES] + [RECORD_QC_VERDICT_TOOL]

QC_VERDICT_QUESTION = """\
Assess whether this dataset's fixed-core QC and clustering results are usable for downstream \
analysis -- the gate that runs BEFORE any known-biology checklist or novel-finding work, so a bad \
run doesn't get expensive follow-up work built on top of it.

First call `get_qc_summary` and read the real numbers (genes/UMIs/mito/doublets per cell, ATAC \
fragments/FRiP/TSS enrichment/nucleosome signal, cluster counts, RNA-ATAC cluster agreement) -- \
don't guess at them. Then spot-check that real biological structure is actually recoverable, not \
just that the QC numbers look fine on paper: call `list_clusters`, pick a couple of clusters and \
check `top_cluster_markers`/`top_gene_activity_markers` on them, and run `cross_modal_marker_check` \
on a marker gene you notice, to confirm at least some of the clustering cross-validates between \
modalities. A moderately elevated doublet rate is a caveat to name, not by itself grounds for \
"do not use," if everything else looks normal -- but an unexpected fragment depth, TSS enrichment, \
cluster agreement, or identity-recovery number IS a real finding to weigh on its own merits.

End by calling `record_qc_verdict` exactly once with your status and specific, numbers-grounded \
reasons.
"""


@dataclass
class QCVerdict:
    status: str
    reasons: list[str] = field(default_factory=list)


def generate_qc_verdict(
    mdata, qc_summary: str, dataset_context: str | None = None, model: str | None = None, max_turns: int = 15,
) -> tuple[QCVerdict, AgentRunResult]:
    """Runs the scoped QC-verdict agent and extracts its `record_qc_verdict`
    call into a typed `QCVerdict`. `dataset_context` defaults to
    `OWN_DATA_CONTEXT` (see `checklist_generator.generate_checklist`'s
    docstring for why that's the right default for arbitrary data).

    If the agent somehow ends without calling `record_qc_verdict` (the
    `require_tool_call` nudge should make this rare, not impossible),
    defaults to "usable with caveats" with an explicit note -- the
    cautious middle ground, not a silent "usable" nor a report-blocking
    "do not use" based on nothing.
    """
    result = run_agent(
        QC_VERDICT_QUESTION, model=model, max_turns=max_turns, mdata=mdata, qc_summary=qc_summary,
        tools=QC_VERDICT_TOOLS, dataset_context=dataset_context if dataset_context is not None else OWN_DATA_CONTEXT,
        require_tool_call="record_qc_verdict",
    )
    calls = [tc for tc in result.tool_calls if tc["name"] == "record_qc_verdict" and not tc.get("is_error", False)]
    if calls:
        verdict_input = calls[-1]["input"]
        verdict = QCVerdict(
            status=verdict_input.get("status", "usable with caveats"),
            reasons=list(verdict_input.get("reasons", [])),
        )
    else:
        logger.warning("generate_qc_verdict: agent never called record_qc_verdict, defaulting to caution")
        verdict = QCVerdict(
            status="usable with caveats",
            reasons=["No verdict was recorded by the agent; defaulting to this cautious status rather than assuming the data is clean."],
        )
    logger.info("QC verdict: %s (%d reason(s)), cost=$%.5f", verdict.status, len(verdict.reasons), result.estimated_cost_usd)
    return verdict, result
