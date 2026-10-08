"""Orchestrates the full 5-section comprehensive report per dataset:

1. Dataset & fixed-core summary
2. Gene activity + chromVAR + cross-modal cell-type-call validation
3. Known-biology checklist via literature RAG (identity discovered, not told)
4. Fault-injection eval
5. Novel findings + adversarial judging + limitations

This module is orchestration and markdown rendering ONLY -- every actual
analysis step reuses an already-built, already-validated piece
(loader_selector, checklist_generator, novelty, fault-injection/scoring,
cross_modal_validation). No new analysis logic lives here.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from multiome_agent.agent.checklist_generator import ChecklistItem, generate_checklist
from multiome_agent.agent.novelty import JudgedFinding, check_novelty_negative_control, judge_novel_findings, propose_novel_findings
from multiome_agent.agent.prompts import PBMC_DATASET_CONTEXT, OWN_DATA_CONTEXT
from multiome_agent.config import REPO_ROOT
from multiome_agent.core.clustering import RNA_CLUSTER_KEY
from multiome_agent.core.cross_modal_validation import systematic_cross_modal_sweep
from multiome_agent.data.dispatch import load_fixed_core_via_agent_decision
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

# Both current datasets are control-only, verified by `check_for_condition_groups`
# during real runs of this pipeline -- no dataset actually has a drug/
# condition axis yet, so there is nothing to render in that branch today.
# It stays available to the agent (see agent/loop.py's tool list) for
# whenever a dataset that does have one is pointed at this pipeline.
DRUG_BRANCH_NOTE = (
    "This dataset has no drug/treatment condition axis (checked by column "
    "name via `check_for_condition_groups`, not assumed) -- every cell is a "
    "control. The pipeline supports condition-level analysis (per-arm QC via "
    "`condition_group_qc`, condition-specific literature RAG) for a future "
    "dataset that does have one; it's simply not exercised here."
)


@dataclass
class ReportCost:
    steps: dict = field(default_factory=dict)

    def add(self, step: str, cost: float) -> None:
        self.steps[step] = self.steps.get(step, 0.0) + cost

    @property
    def total(self) -> float:
        return sum(self.steps.values())

    def render(self) -> str:
        lines = [f"- {step}: ${cost:.4f}" for step, cost in self.steps.items()]
        lines.append(f"- **Total: ${self.total:.4f}**")
        return "\n".join(lines)


def _dataset_summary_text(source: str, mdata) -> str:
    if source == "tenx-cell-ranger":
        from multiome_agent.agent.qc_summary import fixed_core_summary, format_qc_summary
        return format_qc_summary(fixed_core_summary(mdata))
    from multiome_agent.agent.shareseq_qc_summary import format_shareseq_qc_summary, shareseq_fixed_core_summary
    return format_shareseq_qc_summary(shareseq_fixed_core_summary(mdata))


def _render_section1(source: str, dataset_summary_text: str, loader_decision) -> str:
    label = "public 10x PBMC multiome" if source == "tenx-cell-ranger" else "private multi-cell-line multiome"
    return (
        "## 1. Dataset & fixed-core summary\n\n"
        f"**Data source:** {label} data. Loading strategy was determined by the agent itself "
        f"from real file structure, not pre-specified: it decided `{loader_decision.strategy}` "
        f"-- \"{loader_decision.reasoning}\"\n\n"
        f"{dataset_summary_text}\n"
    )


def _render_section2(cross_modal_checks: list[dict], n_rna_clusters: int) -> str:
    lines = [
        "## 2. Gene activity, chromVAR motif deviations, and cross-modal validation\n\n",
        "Gene activity scores and chromVAR-style motif deviations are part of the fixed-core "
        "pipeline (computed once, cached, reused here). The table below is a systematic sweep "
        "(zero extra LLM cost -- these are deterministic Python tool wrappers, run directly, "
        "not routed through an agent turn) of the textbook ArchR/Signac cross-modal cell-type-"
        "call validation over EVERY RNA cluster: its own top-scoring marker gene is checked "
        "against its real cross-modal ATAC partner (matched by cell overlap, not by clusters "
        "coincidentally sharing an integer label). A cluster whose top-scoring gene doesn't "
        "clear the standard significance threshold (padj<0.05, logFC>1.0) is skipped, not "
        "forced with a non-significant marker.\n\n",
        "**How \"Matched ATAC cluster\" is determined:** since RNA and ATAC come from the same "
        "cells, every cell has both an RNA cluster label and an independent ATAC cluster label "
        "(the two clusterings are computed separately, so cluster numbers between them don't "
        "inherently mean anything -- \"RNA cluster 4\" and \"ATAC cluster 4\" aren't related just "
        "because they share a number). To find the real cross-modal partner of an RNA cluster:\n"
        "1. Take every cell in that RNA cluster.\n"
        "2. Look up which ATAC cluster each of those same cells landed in.\n"
        "3. Build a contingency table (RNA cluster x ATAC cluster cell counts) and take the ATAC "
        "cluster with the most overlapping cells -- i.e., the mode.\n\n",
        "| Gene | RNA marker of cluster | Matched ATAC cluster | Gene-activity confirms |\n",
        "|---|---|---|---|\n",
    ]
    if not cross_modal_checks:
        lines.append("| *(no cluster had a significant top marker this run)* | | | |\n")
    for c in cross_modal_checks:
        lines.append(
            f"| {c.get('gene', '?')} | {c.get('rna_cluster', '—')} | "
            f"{c.get('matched_atac_cluster', '—')} | {c.get('gene_activity_confirms', '—')} |\n"
        )
    n_confirmed = sum(1 for c in cross_modal_checks if c.get("gene_activity_confirms") is True)
    n_skipped = n_rna_clusters - len(cross_modal_checks)
    lines.append(
        f"\n{n_confirmed}/{len(cross_modal_checks)} clusters' top markers cross-validate between "
        f"modalities ({n_skipped}/{n_rna_clusters} RNA cluster(s) skipped: no marker cleared the "
        "significance threshold). Per CLAUDE.md's design principle, discordance here is expected, "
        "not a bug -- gene activity is a noisier, indirect accessibility proxy, and RNA-ATAC "
        "agreement is not assumed to be perfect.\n"
    )
    return "".join(lines)


_IDENTITY_LINE_RE = re.compile(r"dataset identity\**\s*:\**\s*(.+)", re.IGNORECASE)


def _truncate_at_boundary(text: str, max_chars: int) -> str:
    """Plain `text[:max_chars]` can (and once did, on a real shareseq-multi-cell-lines
    run whose identity blurb included an 8-row markdown table) slice
    through the middle of a table row, shipping broken markdown in the
    committed report. Prefer cutting at a paragraph break, then a line
    break, and only fall back to a raw character cut if neither exists
    before the limit.
    """
    if len(text) <= max_chars:
        return text
    cut = text.rfind("\n\n", 0, max_chars)
    if cut == -1:
        cut = text.rfind("\n", 0, max_chars)
    if cut == -1:
        cut = max_chars
    return text[:cut].rstrip() + "..."


def _extract_identity_blurb(answer: str, max_chars: int = 900) -> str:
    """The checklist-generation agent's final answer commonly restates the
    whole checklist as a structured breakdown (often after a leading "---"
    and a "## Summary" heading) -- that table is redundant with (and less
    exact than) the `items` list rendered right below this blurb. Verified
    against two real runs' actual final-answer text, not a guessed shape:
    both led with "---\n\n## Summary\n\n**Dataset Identity**: <statement>."
    Prefer pulling that explicit line; fall back to "prose before the first
    heading" for an answer that states identity directly with no heading at
    all (also seen in a real run); fall back to a pointer if neither yields
    anything -- never dump the full raw multi-paragraph answer.
    """
    text = answer.strip()

    m = _IDENTITY_LINE_RE.search(text)
    if m:
        return _truncate_at_boundary(m.group(1).strip(), max_chars)

    text = re.sub(r"^\s*-{3,}\s*\n+", "", text)  # drop a leading horizontal-rule line
    text = re.sub(r"^\s*#{1,6}\s.*\n+", "", text)  # drop a leading heading line too (e.g. "## Summary")
    cut = re.search(r"\n\s*#{1,6}\s", text)
    if cut:
        text = text[: cut.start()].strip()
    if not text or re.fullmatch(r"-{3,}", text):
        return "(see this run's reasoning trail for the agent's full identity determination)"
    return _truncate_at_boundary(text, max_chars)


def _render_section3(items: list[ChecklistItem], identity_answer: str) -> str:
    lines = [
        "## 3. Known-biology checklist (literature RAG, identity determined by the agent)\n\n",
        f"**Identity determination:** {_extract_identity_blurb(identity_answer)}\n\n",
        "Each item below required the agent to: search PubMed, fetch and actually read a real "
        "abstract (not just a title), and verify with a real tool call that the gene/motif is "
        "usable in this dataset before recording it.\n\n",
    ]
    labels = {
        "rna_marker": "RNA markers", "motif": "Motifs",
        "peak_to_gene": "Distal peak-to-gene links", "regulon_target": "TF-target-gene regulon links",
    }
    confirmed = [it for it in items if it.confirmed_present_in_data]
    rejected = [it for it in items if not it.confirmed_present_in_data]

    lines.append("## Confirmed by data\n\n")
    by_category = {"rna_marker": [], "motif": [], "peak_to_gene": [], "regulon_target": []}
    for item in confirmed:
        by_category.setdefault(item.category, []).append(item)
    for cat, cat_items in by_category.items():
        lines.append(f"### {labels.get(cat, cat)}\n\n")
        if not cat_items:
            lines.append("*(no grounded, data-present candidate found this run)*\n\n")
            continue
        for it in cat_items:
            cell_line_note = f" (confirmed in: {', '.join(it.cell_lines)})" if it.cell_lines else ""
            lines.append(f"- **{it.gene_or_motif}**: {it.claim} (PMID {it.pmid}, *{it.journal}*, {it.year}){cell_line_note}\n")
        lines.append("\n")

    lines.append("## Rejected by data\n\n")
    if not rejected:
        lines.append("*(no candidate this run clearly failed its own data-verification check)*\n\n")
    else:
        lines.append(
            "Real literature-backed hypotheses whose own data-verification tool call clearly "
            "contradicted the claim (not significant, wrong-sign, or absent from this dataset) -- "
            "not candidates that simply weren't needed once enough confirmed items were found:\n\n"
        )
        for it in rejected:
            lines.append(
                f"- **{it.gene_or_motif}** ({labels.get(it.category, it.category)}): {it.claim} "
                f"(PMID {it.pmid}, *{it.journal}*, {it.year})\n"
            )
        lines.append("\n")
    return "".join(lines)


FAULT_INJECTION_MODELS = ("claude-haiku-4-5", "claude-sonnet-5", "claude-opus-5", "claude-fable-5")


def _render_section4(results_by_model: dict[str, dict], narrative_model: str) -> str:
    """Renders a cross-model fault-detection comparison table (Haiku/Sonnet/
    Opus/Fable, per the original CLAUDE.md spec for this section
    specifically) plus one model's full per-scenario answers for
    qualitative detail -- showing all 4 models' full text would bloat the
    report; the comparison table is the actual cross-model content.
    """
    lines = [
        "## 4. Fault-injection eval (compared across models)\n\n",
        "| Model | Faults detected | Correct diagnosis | False alarms | Cost/run |\n",
        "|---|---|---|---|---|\n",
    ]
    for model, results in results_by_model.items():
        scen = results["scenario_results"]
        fault_scen = [r for r in scen if r["fault_type"] != "clean_control"]
        clean_scen = [r for r in scen if r["fault_type"] == "clean_control"]
        detected = sum(1 for r in fault_scen if r["detected"])
        diag = sum(1 for r in fault_scen if r["diagnosis_matches"])
        false_alarms = sum(1 for r in clean_scen if r["false_alarm"])
        avg_cost = results["total_cost_usd"] / len(scen) if scen else 0.0
        lines.append(
            f"| {model} | {detected}/{len(fault_scen)} | {diag}/{len(fault_scen)} | "
            f"{false_alarms}/{len(clean_scen)} | ${avg_cost:.5f} |\n"
        )
    lines.append(
        "\n\"Faults detected\"/\"Correct diagnosis\" use heuristic free-text classifiers on each "
        "model's open-ended answer, not exact ground-truth string matching.\n"
    )

    lines.append(f"\n### Per-scenario answers ({narrative_model})\n")
    for r in results_by_model[narrative_model]["scenario_results"]:
        lines.append(f"\n**{r['fault_type']}** ({r['description']})\n\n{r['answer']}\n")
    return "".join(lines)


def _render_judged_finding(j: JudgedFinding) -> str:
    v = j.verdict
    # `record_judger_verdict`'s "required" fields are a hint to the
    # model, not server-enforced -- a real Opus run once called the
    # tool without a "verdict" key at all (KeyError, crashing the
    # report after every expensive step had already run and been
    # paid for). Use .get() throughout so a partially-filled verdict
    # still renders something honest instead of crashing.
    verdict_str = v.get("verdict", "NO VERDICT RECORDED") if v else "NO VERDICT RECORDED"
    cell_line_note = f" *(tested in: {', '.join(j.cell_lines)})*" if j.cell_lines else ""
    out = f"### Finding (confidence: {j.confidence}){cell_line_note}\n\n{j.finding}\n\n**Evidence:** {j.evidence}\n\n"
    out += f"**Judger verdict: {verdict_str}**\n\n"
    if v:
        out += (
            f"- Likely artifact: {v.get('likely_artifact', '(not provided)')} -- "
            f"{v.get('artifact_reasoning', '(not provided)')}\n"
            f"- Already known: {v.get('already_known', '(not provided)')} -- "
            f"{v.get('novelty_reasoning', '(not provided)')}\n\n"
        )
    return out


def _render_unjudged_concern(f: dict) -> str:
    return (
        f"### Finding, not adversarially judged (confidence: {f.get('confidence', '?')})\n\n"
        f"{f.get('finding', '')}\n\n**Evidence:** {f.get('evidence', '')}\n\n"
        f"*Classified by the model itself as `no_signal_or_concern` -- reports an absence of "
        f"signal or a methodological/tooling observation rather than a positive discovery claim, "
        f"so it wasn't sent to the adversarial Judger (there's no positive claim to stress-test).*\n\n"
    )


def _render_section5(
    judged: list[JudgedFinding], negative_control: dict, limitations: list[str]
) -> str:
    lines = ["## 5. Novel findings, adversarial judging, and limitations\n\n"]
    if not judged:
        lines.append("No novel findings were proposed this run (a legitimate outcome, not a failure).\n\n")
    for j in judged:
        lines.append(_render_judged_finding(j))

    lines.append("### Negative control: novel-finding proposal under shuffled RNA-ATAC pairing\n\n")
    lines.append(
        f"With no real cross-modal relationship left (100% shuffled pairing), the same "
        f"novel-finding-proposal step was run again. It proposed "
        f"{negative_control['n_findings_proposed']} finding(s) "
        f"(0 is the well-behaved outcome; any finding the model itself classifies as a genuine "
        f"positive relationship claim gets adversarially judged below -- surviving judging here "
        f"would indicate real hallucination from noise; findings classified as reporting an "
        f"absence of signal or a tooling concern aren't judged, since there's no positive claim "
        f"to stress-test).\n\n"
    )
    for j in negative_control.get("judged_positive_findings", []):
        lines.append(_render_judged_finding(j))
    for f in negative_control.get("no_signal_or_concern_findings", []):
        lines.append(_render_unjudged_concern(f))

    lines.append("### Limitations\n\n")
    for lim in limitations:
        lines.append(f"- {lim}\n")
    return "".join(lines)


_GENERAL_LIMITATIONS = [
    "\"Faults detected\"/\"diagnosis matches\" use heuristic free-text classifiers on the "
    "agent's open-ended answer, not exact ground-truth string matching.",
    "Known-biology checklist items are generated fresh each run via literature RAG, not scored "
    "against a fixed pre-written ground-truth file -- recall/precision against a fixed checklist "
    "is a different, complementary evaluation this report doesn't repeat.",
    "Gene activity is a noisier, indirect accessibility proxy than direct RNA counts; "
    "cross-modal disagreement on a real marker is expected, not necessarily a data-quality issue.",
    "This report reflects a single run on one model; consistency across repeated runs (or across "
    "models) is a separate axis this report doesn't cover.",
]

_SHARESEQ_ONLY_LIMITATIONS = [
    "ATAC downsampling and doublet injection (used for the public dataset) could not be reused "
    "as-is here: this dataset has no raw fragments file (downsampling needs one) and RNA `.X` "
    "has no raw counts layer (doublet injection needs one) -- a real format constraint, not a "
    "design choice. The fault-injection section uses cell-line-label-swap, mixed-samples, and "
    "shuffled-RNA-ATAC-pairing instead.",
]


def generate_report(source: str, model: str, out_path: Path | None = None) -> dict:
    """Generates the full 5-section report for `source` ('tenx-cell-ranger' or
    'shareseq-multi-cell-lines') on `model`, writes it to `out_path` (default:
    reports/examples/{source}-report.md), and returns {report_text, cost, ...}
    for programmatic inspection.
    """
    from multiome_agent.core.condition_detection import detect_condition_groups

    cost = ReportCost()
    dataset_context = PBMC_DATASET_CONTEXT if source == "tenx-cell-ranger" else OWN_DATA_CONTEXT

    mdata, loader_decision = load_fixed_core_via_agent_decision(source, model=model)
    cost.add("loader_decision", loader_decision.cost_usd)
    dataset_summary_text = _dataset_summary_text(source, mdata)

    condition_check = detect_condition_groups(mdata)

    # Section 3 mechanism: identity discovery + literature-RAG checklist.
    items, checklist_result = generate_checklist(mdata, dataset_summary_text, dataset_context, model=model)
    cost.add("checklist_generation", checklist_result.estimated_cost_usd)

    checklist_summary = "\n".join(f"- ({i.category}) {i.claim} [PMID {i.pmid}]" for i in items) or "(no checklist items)"

    # Section 5 mechanism: propose + adversarially judge.
    findings, novelty_result = propose_novel_findings(mdata, dataset_summary_text, dataset_context, checklist_summary, model=model)
    cost.add("novelty_proposal", novelty_result.estimated_cost_usd)
    judged = judge_novel_findings(findings, mdata, dataset_summary_text, dataset_context, model=model)
    cost.add("judging", sum(j.judge_cost_usd for j in judged))

    # Section 4: fault-injection eval, compared across all 4 Claude models
    # (Haiku/Sonnet/Opus/Fable) -- scenarios (including the one real ~17-min
    # snapatac2 downsampling rerun for the tenx-cell-ranger dataset) are built
    # ONCE and reused across every model, not rebuilt per model.
    if source == "tenx-cell-ranger":
        from multiome_agent.eval.eval_harness import build_scenarios, run_model_eval
        scenarios = build_scenarios(mdata)
        fault_injection_results_by_model = {m: run_model_eval(m, scenarios, mdata) for m in FAULT_INJECTION_MODELS}
    else:
        from multiome_agent.eval.shareseq_eval_harness import build_shareseq_scenarios, run_shareseq_model_eval
        scenarios = build_shareseq_scenarios(mdata)
        fault_injection_results_by_model = {m: run_shareseq_model_eval(m, scenarios, mdata) for m in FAULT_INJECTION_MODELS}
    cost.add("fault_injection", sum(r["total_cost_usd"] for r in fault_injection_results_by_model.values()))

    # Negative control: propose novel findings against an ACTUALLY shuffled
    # (not just a doctored qc_summary) RNA-ATAC pairing, recomputed fresh --
    # every real tool call the agent makes (cross_modal_marker_check,
    # top_cluster_markers) must see the true faulted state, not the clean one.
    if source == "tenx-cell-ranger":
        from multiome_agent.eval.fault_injection import shuffle_rna_atac_pairing
        from multiome_agent.agent.qc_summary import fixed_core_summary, format_qc_summary
        shuffled_mdata, _ = shuffle_rna_atac_pairing(mdata, fraction=1.0, seed=0)
        shuffled_qc = format_qc_summary(fixed_core_summary(shuffled_mdata))
    else:
        from multiome_agent.eval.shareseq_fault_injection import shuffle_rna_atac_pairing
        from multiome_agent.agent.shareseq_qc_summary import format_shareseq_qc_summary, shareseq_fixed_core_summary
        shuffled_mdata, _ = shuffle_rna_atac_pairing(mdata, fraction=1.0, seed=0)
        shuffled_qc = format_shareseq_qc_summary(shareseq_fixed_core_summary(shuffled_mdata))
    negative_control = check_novelty_negative_control(shuffled_mdata, shuffled_qc, dataset_context, model=model)
    cost.add("negative_control", negative_control["cost_usd"])

    cross_modal_checks = systematic_cross_modal_sweep(mdata)
    n_rna_clusters = mdata.mod["rna"].obs[RNA_CLUSTER_KEY].nunique()

    limitations = list(_GENERAL_LIMITATIONS)
    if source == "shareseq-multi-cell-lines":
        limitations.extend(_SHARESEQ_ONLY_LIMITATIONS)
    if condition_check["is_control_only"]:
        limitations.append(DRUG_BRANCH_NOTE)

    sections = [
        _render_section1(source, dataset_summary_text, loader_decision),
        _render_section2(cross_modal_checks, n_rna_clusters),
        _render_section3(items, checklist_result.answer),
        _render_section4(fault_injection_results_by_model, model),
        _render_section5(judged, negative_control, limitations),
    ]
    report_text = (
        f"# Multiome QC & Hypothesis Agent report ({source})\n\n"
        f"*Model: `{model}` · total cost: ${cost.total:.4f}*\n\n"
        + "\n---\n\n".join(sections)
        + f"\n\n---\n\n## Cost breakdown\n\n{cost.render()}\n"
    )

    out_path = out_path or (REPO_ROOT / "reports" / "examples" / f"{source}-report.md")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report_text)
    logger.info("Wrote %s report to %s (total cost $%.4f)", source, out_path, cost.total)

    return {"report_text": report_text, "cost": cost, "out_path": str(out_path)}
