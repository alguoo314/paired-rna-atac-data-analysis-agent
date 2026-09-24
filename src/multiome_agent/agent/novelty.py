"""Novel-finding proposal + adversarial judging.

Two agent roles, same structured-tool-call pattern as checklist_generator.py:
1. Propose up to 3 novel findings -- specific, checkable patterns in the
   data that go beyond the already-confirmed known-biology checklist.
2. Judge each one independently and adversarially: is it likely an
   artifact? Does further literature RAG show it's already known? A
   finding surviving both checks is a stronger claim than one that never
   faced either. Per spec, it's fine -- expected, even -- for every
   proposed finding to be struck down.

Also provides the negative-control check the shareseq-multi-cell-lines-data fault-injection
section needs: run the SAME proposal step against a shuffled-RNA-ATAC-pairing
(no real cross-modal relationship left) state and confirm it doesn't
hallucinate a spurious cross-modal finding.
"""

from __future__ import annotations

from dataclasses import dataclass

from multiome_agent.agent.loop import TOOLS, AgentRunResult, run_agent
from multiome_agent.agent.prompts import OWN_DATA_CONTEXT
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

RECORD_NOVEL_FINDING_TOOL = {
    "name": "record_novel_finding",
    "description": (
        "Record ONE candidate novel finding that survived YOUR OWN literature check: a specific, "
        "checkable pattern in this data that goes beyond the already-confirmed known-biology "
        "checklist, AND that a real search_pubmed (+ fetch_pubmed_abstracts on any promising hit) "
        "you just ran did NOT already find clearly reported. If your literature check turns up "
        "that exact relationship, do NOT call this -- silently move on to a different candidate "
        "instead. You may try up to 6 distinct candidates this way, but stop calling this tool "
        "once you've recorded 3 that survived their own check."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "finding": {"type": "string", "description": "The specific claim, e.g. 'motif X accessibility correlates with gene Y expression specifically in cluster Z'"},
            "evidence": {"type": "string", "description": "The specific real tool result(s) supporting this, with actual numbers"},
            "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
            "finding_type": {
                "type": "string",
                "enum": ["positive_relationship", "no_signal_or_concern"],
                "description": (
                    "'positive_relationship': this claims a specific, checkable pattern or "
                    "relationship IS present in the data -- a real discovery claim, which gets "
                    "adversarially judged. 'no_signal_or_concern': this instead reports that an "
                    "expected relationship is ABSENT (e.g. correctly finding nothing under a "
                    "negative-control/shuffled state) or raises a methodological/tooling concern "
                    "without itself claiming novel positive biology (e.g. a limitation in how a "
                    "tool matches clusters) -- these are NOT adversarially judged, since there's "
                    "no positive claim to stress-test for hallucination."
                ),
            },
        },
        "required": ["finding", "evidence", "confidence", "finding_type"],
    },
}
NOVELTY_TOOLS = TOOLS + [RECORD_NOVEL_FINDING_TOOL]

RECORD_JUDGER_VERDICT_TOOL = {
    "name": "record_judger_verdict",
    "description": "Record your final verdict on the one proposed novel finding you were asked to judge. Call exactly once.",
    "input_schema": {
        "type": "object",
        "properties": {
            "likely_artifact": {"type": "boolean", "description": "True if this is likely a technical/statistical artifact rather than real biology"},
            "artifact_reasoning": {"type": "string"},
            "already_known": {"type": "boolean", "description": "True if your own literature search found this already well-established"},
            "novelty_reasoning": {"type": "string"},
            "verdict": {"type": "string", "enum": ["survives", "struck_down"]},
        },
        "required": ["likely_artifact", "artifact_reasoning", "already_known", "novelty_reasoning", "verdict"],
    },
}
JUDGER_TOOLS = TOOLS + [RECORD_JUDGER_VERDICT_TOOL]


def _novelty_question(checklist_summary: str, is_negative_control: bool = False) -> str:
    covered = checklist_summary or "(no checklist items available)"
    classify_note = (
        """\
Classify each finding you record with `finding_type`: "positive_relationship" if you're claiming \
a specific pattern IS present (this will be adversarially judged), or "no_signal_or_concern" if \
you're instead reporting that an expected relationship is ABSENT or flagging a methodological/\
tooling concern rather than claiming discovered biology. Given the RNA-ATAC pairing here is \
COMPLETELY SHUFFLED (a negative control -- no real cross-modal relationship can exist), the \
honest, well-behaved answer for most candidates is "no_signal_or_concern" or no finding at all; \
only use "positive_relationship" if you genuinely believe you've found a specific pattern that \
would survive scrutiny even knowing the pairing is fake, and be able to defend that against a \
skeptical judge.
""" if is_negative_control else (
        """\
Classify each finding you record with `finding_type`: "positive_relationship" for a genuine \
discovery claim (the normal case here), or "no_signal_or_concern" only if you're instead \
reporting an absence of an expected relationship or a methodological/tooling concern rather than \
claiming new biology.
"""
    ))
    return f"""\
Based on everything you've learned about this dataset's real identity, and these \
already-established literature findings you (or a prior step) confirmed:

{covered}

Try up to 6 candidate NOVEL findings: specific, checkable patterns in THIS data that go beyond \
what's already covered above -- not a restatement of a known-biology checklist item. For EACH \
candidate, before deciding whether to record it: (1) gather the real tool evidence for it (actual \
numbers, not guessed), then (2) run at least one `search_pubmed` call -- and `fetch_pubmed_abstracts` \
on any promising hit -- to check whether this SAME specific pattern is already reported. If your \
own literature check already shows it's well-established, SELF-REJECT that candidate silently: do \
NOT call `record_novel_finding` for it, just move on and try a different candidate. Only call \
`record_novel_finding` for a candidate that survives its own literature check, with an honest \
confidence level (high/medium/low). Stop once you've recorded 3 findings that survived this way, \
or once you've tried 6 distinct candidates, whichever comes first. It is completely fine -- \
expected, even -- to end up with fewer than 3 recorded findings (even zero) if every candidate you \
tried turned out to already be known -- don't force a weaker or already-known candidate just to \
fill the quota.

If this dataset pools multiple distinct cell lines/lineages, don't concentrate all your candidates \
on the same single line -- try to span them across different lines chosen at random (one about \
line A, another about line B, another about line C), each still as specific and well-evidenced as \
if you'd focused on one line the whole time.

{classify_note}"""


def _judger_question(finding: str, evidence: str) -> str:
    return f"""\
You are acting as a skeptical, adversarial judge of ONE proposed "novel finding" from a multiome \
analysis. Your default assumption is that it is EITHER an artifact OR already known in the \
literature -- only let it survive if you genuinely can't find a reason to reject it after real \
investigation.

Proposed finding: {finding}
Evidence originally cited for it: {evidence}

Investigate independently, using your own tool calls (don't just take the original evidence at \
face value):
1. Artifact check: look for alternative, boring explanations using the data-analysis tools \
available to you -- e.g. a small sample size in the relevant cluster, a correlation plausibly \
driven by a confound, or an effect that doesn't independently replicate when checked a different \
way (e.g. via cross_modal_marker_check).
2. Novelty check: use search_pubmed then fetch_pubmed_abstracts (actually read the abstract, \
don't just check titles) to look for existing literature that already reports this same \
relationship.

Then call `record_judger_verdict` exactly once with your verdict.
"""


@dataclass
class JudgedFinding:
    finding: str
    evidence: str
    confidence: str
    verdict: dict | None
    judge_cost_usd: float
    finding_type: str = "positive_relationship"


def propose_novel_findings(
    mdata, qc_summary: str, dataset_context: str | None = None, checklist_summary: str = "",
    model: str | None = None, max_turns: int = 45, is_negative_control: bool = False,
) -> tuple[list[dict], AgentRunResult]:
    """`dataset_context` defaults to `OWN_DATA_CONTEXT` (see
    `checklist_generator.generate_checklist`'s docstring for why). Default
    `max_turns` raised from the original 25: trying up to 6 candidates, each
    needing its own tool-evidence gathering AND a literature self-check
    before the model decides whether to record it, needs real headroom.
    """
    result = run_agent(
        _novelty_question(checklist_summary, is_negative_control=is_negative_control), model=model,
        max_turns=max_turns, mdata=mdata, qc_summary=qc_summary, tools=NOVELTY_TOOLS,
        dataset_context=dataset_context if dataset_context is not None else OWN_DATA_CONTEXT,
    )
    findings = [tc["input"] for tc in result.tool_calls if tc["name"] == "record_novel_finding" and not tc.get("is_error", False)]
    for f in findings:
        f.setdefault("finding_type", "positive_relationship")  # schema "required" isn't server-enforced
    # Deterministic cap, not relied on the model's own discipline (same
    # pattern as checklist_generator's dedup): the prompt asks the model to
    # stop recording once it has 3 findings that survived their own
    # literature self-check, but don't trust that alone. Only caps
    # "positive_relationship" findings -- "no_signal_or_concern" findings
    # aren't judged downstream anyway, so there's no reason to cap those.
    capped, n_positive = [], 0
    for f in findings:
        if f.get("finding_type") == "positive_relationship":
            n_positive += 1
            if n_positive > 3:
                continue
        capped.append(f)
    if n_positive > 3:
        logger.warning("Model recorded %d positive_relationship findings; keeping only the first 3", n_positive)
    findings = capped
    logger.info("Proposed %d novel finding(s), cost=$%.5f", len(findings), result.estimated_cost_usd)
    return findings, result


def judge_novel_findings(
    findings: list[dict], mdata, qc_summary: str, dataset_context: str | None = None,
    model: str | None = None, max_turns: int = 20,
) -> list[JudgedFinding]:
    """`dataset_context` defaults to `OWN_DATA_CONTEXT` (see
    `checklist_generator.generate_checklist`'s docstring for why)."""
    dataset_context = dataset_context if dataset_context is not None else OWN_DATA_CONTEXT
    judged = []
    for f in findings:
        result = run_agent(
            _judger_question(f["finding"], f["evidence"]), model=model, max_turns=max_turns, mdata=mdata,
            qc_summary=qc_summary, tools=JUDGER_TOOLS, dataset_context=dataset_context,
            require_tool_call="record_judger_verdict",
        )
        verdict_calls = [tc for tc in result.tool_calls if tc["name"] == "record_judger_verdict" and not tc.get("is_error", False)]
        verdict = verdict_calls[-1]["input"] if verdict_calls else None
        # `record_judger_verdict`'s "verdict" enum is schema-"required" but not
        # server-enforced -- a real Opus run called the tool with thorough,
        # correct artifact_reasoning/already_known but simply omitted the
        # "verdict" field itself. Rather than pay for another nudge round
        # that risks the same omission again, backfill deterministically from
        # the two booleans the model DID commit to: the judging prompt's own
        # rule is "only survives if you can't find a reason to reject it",
        # i.e. struck_down iff likely_artifact or already_known.
        if verdict is not None and "verdict" not in verdict:
            verdict["verdict"] = "struck_down" if verdict.get("likely_artifact") or verdict.get("already_known") else "survives"
        judged.append(JudgedFinding(
            finding=f["finding"], evidence=f["evidence"], confidence=f["confidence"],
            verdict=verdict, judge_cost_usd=result.estimated_cost_usd,
            finding_type=f.get("finding_type", "positive_relationship"),
        ))
        logger.info("Judged finding %r -> %s", f["finding"][:60], verdict.get("verdict") if verdict else "NO VERDICT")
    return judged


def check_novelty_negative_control(
    shuffled_mdata, shuffled_qc_summary: str, dataset_context: str | None = None, model: str | None = None
) -> dict:
    """Runs `propose_novel_findings` against a shuffled-RNA-ATAC-pairing
    state (no real cross-modal relationship left) -- a well-behaved
    analysis should propose zero cross-modal findings, or none at all,
    rather than hallucinating a spurious relationship from noise.

    Only findings the model itself classified as `finding_type ==
    "positive_relationship"` (a genuine discovery claim, the failure mode
    this check exists to catch) go through the adversarial Judger --
    "no_signal_or_concern" findings (correctly reporting an absence of
    signal, or a methodological/tooling observation) aren't positive
    claims to stress-test, so judging them would just spend money without
    testing anything. Real example that motivated this split: a genuine
    negative-control run once proposed a finding that the cross-modal
    cluster-matching TOOL itself is degenerate (a real, shuffle-independent
    observation about the pipeline, not a hallucinated relationship) -- no
    reason to burn a Judger call on that.
    """
    findings, result = propose_novel_findings(
        shuffled_mdata, shuffled_qc_summary, dataset_context, model=model, is_negative_control=True,
    )
    positive = [f for f in findings if f.get("finding_type") == "positive_relationship"]
    concerns = [f for f in findings if f.get("finding_type") != "positive_relationship"]
    judged = judge_novel_findings(positive, shuffled_mdata, shuffled_qc_summary, dataset_context, model=model) if positive else []
    judge_cost = sum(j.judge_cost_usd for j in judged)
    logger.info(
        "Negative control: %d finding(s) (%d positive_relationship -> judged, %d no_signal_or_concern -> not judged)",
        len(findings), len(positive), len(concerns),
    )
    return {
        "n_findings_proposed": len(findings), "findings": findings,
        "judged_positive_findings": judged, "no_signal_or_concern_findings": concerns,
        "cost_usd": result.estimated_cost_usd + judge_cost, "answer": result.answer,
    }
