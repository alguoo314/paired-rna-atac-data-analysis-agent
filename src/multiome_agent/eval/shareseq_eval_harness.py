"""Shareseq-multi-cell-lines-data analog of `eval_harness.py`: ties the shareseq-multi-cell-lines-data fault
types to a real agent quality judgment, scored against ground truth.
Reused directly by the comprehensive report (`agent/report_generator.py`)
§4; the multi-model comparison this file originally also built is retired
(see CLAUDE.md's Evaluation §8 and `run_shareseq_model_eval`'s docstring).

Design choice, matching the tenx-cell-ranger harness exactly and deliberately
avoiding the mistake documented in PROGRESS.md's step 10 (an earlier
version of the tenx-cell-ranger harness handed the model a "baseline vs. current"
comparison directly, inflating detection rates for free): every scenario's
`run_agent` call gets the CLEAN cached `mdata` (neither shareseq-multi-cell-lines fault type
touches RNA/ATAC feature values, only `.obs["cell_line_name"]`, so
`tf_motif_correlation`/`enrich_gene_set` behave identically either way --
reusing the clean object avoids copying the ~1.9GB MuData per scenario too)
but a scenario-specific `qc_summary` string. Every field in that summary is
computed FRESH from that scenario's own (possibly faulted) state -- never
alongside the clean value for comparison, and computed identically across
every faulted scenario (same fields present every time) so the mere
presence/absence of a field can't itself leak which scenario is faulted.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from sklearn.metrics import adjusted_rand_score

from multiome_agent.agent.loop import run_agent
from multiome_agent.agent.shareseq_qc_summary import format_shareseq_qc_summary, shareseq_fixed_core_summary
from multiome_agent.agent.prompts import OWN_DATA_CONTEXT
from multiome_agent.core.clustering import ATAC_CLUSTER_KEY, RNA_CLUSTER_KEY
from multiome_agent.eval.eval_harness import _NEGATION_CUES, _NEGATION_WINDOW, PROBLEM_WORDS, _classify_detected
from multiome_agent.eval.grounding import check_grounding
from multiome_agent.eval.shareseq_fault_injection import shuffle_rna_atac_pairing, swap_cell_line_labels
from multiome_agent.eval.shareseq_scoring import score_shuffled_pairing
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

QUESTION = (
    "This is a processed, pooled multi-cell-line single-cell multiome dataset. Assess its "
    "quality using the QC summary tool. Does anything look wrong? If so, describe the specific "
    "problem and how confident you are; if not, say the data looks clean."
)

CELL_LINE_SWAP_FRACTION = 0.30  # a clear, validated mid-severity signal per
# PROGRESS.md's real degradation curve (0.767 clean -> 0.623 @10% -> 0.358
# @30% -> 0.157 @50%) -- strong enough to be a real, checkable fault without
# being the most extreme (nearly-total) severity level.

SHARESEQ_FAULT_KEYWORDS = {
    "cell_line_label_swap": [
        "cell line", "cell-line", "identity", "genotype", "label", "demultiplex",
        "recovery", "mislabel", "misassign", "swap",
    ],
    "shuffled_rna_atac_pairing": [
        "pairing", "mismatch", "cluster agreement", "cross-modal", " ari", "shuffl", "misalign", "correspond",
    ],
}


def _classify_diagnosis_matches(answer_text: str, fault_type: str) -> bool | None:
    if fault_type == "clean_control":
        return None
    return any(kw in answer_text.lower() for kw in SHARESEQ_FAULT_KEYWORDS[fault_type])


@dataclass
class ShareseqScenarioResult:
    fault_type: str
    severity: float
    description: str
    answer: str
    detected: bool
    diagnosis_matches: bool | None
    false_alarm: bool
    cost: float
    turns: int
    grounding: dict
    reasoning_trail_path: str
    refused: bool = False
    refusal_detail: str | None = None


def _cell_line_recovery(rna_labels, atac_labels, cell_line_labels) -> dict:
    rna_ari = adjusted_rand_score(rna_labels, cell_line_labels.reindex(rna_labels.index))
    atac_ari = adjusted_rand_score(atac_labels, cell_line_labels.reindex(atac_labels.index))
    return {"cell_line_recovery_ari_rna": rna_ari, "cell_line_recovery_ari_atac": atac_ari}


def _format_scenario_summary(
    base_summary: dict, recovery: dict | None = None, cross_modal_ari: float | None = None
) -> str:
    summary = dict(base_summary)
    if recovery is not None:
        summary["cell_line_recovery_ari_rna"] = recovery["cell_line_recovery_ari_rna"]
        summary["cell_line_recovery_ari_atac"] = recovery["cell_line_recovery_ari_atac"]
    if cross_modal_ari is not None:
        summary["cross_modal_ari"] = cross_modal_ari
    return format_shareseq_qc_summary(summary)


def build_shareseq_scenarios(clean_mdata) -> list[dict]:
    """Returns [{fault_type, severity, description, qc_summary}]. Reuses the
    SAME clean `mdata` object for every scenario's `run_agent(mdata=...)`
    call (see module docstring) -- only `qc_summary` varies.
    """
    scenarios = []

    rna = clean_mdata.mod["rna"]
    atac = clean_mdata.mod["atac"]
    clean_summary = shareseq_fixed_core_summary(clean_mdata)
    scenarios.append({
        "fault_type": "clean_control", "severity": 0.0,
        "description": "Unmodified private multi-cell-line data.",
        "qc_summary": _format_scenario_summary(clean_summary),
    })

    swapped_mdata, swap_meta = swap_cell_line_labels(clean_mdata, fraction=CELL_LINE_SWAP_FRACTION, seed=0)
    swapped_labels = swapped_mdata.mod["rna"].obs["cell_line_name"]
    swap_recovery = _cell_line_recovery(rna.obs[RNA_CLUSTER_KEY], atac.obs[ATAC_CLUSTER_KEY], swapped_labels)
    scenarios.append({
        "fault_type": "cell_line_label_swap", "severity": CELL_LINE_SWAP_FRACTION,
        "description": swap_meta.description,
        "qc_summary": _format_scenario_summary(clean_summary, swap_recovery),
    })

    shuffled_mdata, shuffle_meta = shuffle_rna_atac_pairing(clean_mdata, fraction=0.5, seed=0)
    shuffle_score = score_shuffled_pairing(clean_mdata, shuffled_mdata)
    # Cell-line recovery is untouched by this fault (see
    # shuffle_rna_atac_pairing's docstring) -- only cross-modal ARI moves.
    scenarios.append({
        "fault_type": "shuffled_rna_atac_pairing", "severity": 0.5,
        "description": shuffle_meta.description,
        "qc_summary": _format_scenario_summary(clean_summary, cross_modal_ari=shuffle_score["ari_faulted"]),
    })

    return scenarios


def _run_shareseq_scenario(clean_mdata, scenario: dict, model: str | None = None) -> ShareseqScenarioResult:
    result = run_agent(
        QUESTION, model=model, mdata=clean_mdata, qc_summary=scenario["qc_summary"],
        dataset_context=OWN_DATA_CONTEXT,
    )
    grounding = check_grounding(result)
    detected = _classify_detected(result.answer)
    is_clean = scenario["fault_type"] == "clean_control"
    return ShareseqScenarioResult(
        fault_type=scenario["fault_type"], severity=scenario["severity"], description=scenario["description"],
        answer=result.answer, detected=detected,
        diagnosis_matches=_classify_diagnosis_matches(result.answer, scenario["fault_type"]),
        false_alarm=(is_clean and detected), cost=result.estimated_cost_usd, turns=result.turn_count,
        grounding=grounding, reasoning_trail_path=result.reasoning_trail_path,
        refused=result.refused, refusal_detail=result.refusal_detail,
    )


def run_shareseq_model_eval(model: str, scenarios: list[dict], clean_mdata) -> dict:
    """Single-model, no-consistency-repeat runner -- kept as a lighter-weight
    entry point than the comprehensive report's own per-scenario calls.
    The multi-model comparison this file used to also build
    (`run_shareseq_multi_model_comparison`/`build_shareseq_comparison_table`,
    writing `shareseq-multi-cell-lines-model-comparison.md`) is retired -- see CLAUDE.md's
    Evaluation §8; the comprehensive report (`agent/report_generator.py`)
    is the deliverable now, generated once on one model.
    """
    results = [_run_shareseq_scenario(clean_mdata, s, model=model) for s in scenarios]
    return {
        "model": model,
        "scenario_results": [asdict(r) for r in results],
        "total_cost_usd": sum(r.cost for r in results),
    }
