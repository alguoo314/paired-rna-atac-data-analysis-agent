"""Split-half replication + negative-control checks for the shareseq-multi-cell-lines
data (CLAUDE.md's "Novel findings" eval item: "Negative
control: shuffled pairing should yield NO novel relationships. Split-half
replication across random cell halves."). Presupposes the agent can propose
an open-ended "notable pattern" claim at all -- existing tools are narrow
(a named gene's TF-motif check, a given gene list's enrichment), so this
uses a dedicated open-ended question against the QC-summary tool instead.

Reports statistics only, no cell-line/library identity strings, per
CLAUDE.md's data-handling rules for this dataset.
"""

from __future__ import annotations

import numpy as np
from mudata import MuData

from multiome_agent.agent.loop import run_agent
from multiome_agent.agent.shareseq_qc_summary import format_shareseq_qc_summary, shareseq_fixed_core_summary
from multiome_agent.agent.prompts import SHARESEQ_DATASET_CONTEXT
from multiome_agent.core.shareseq_pipeline import run_shareseq_fixed_core
from multiome_agent.data.shareseq_loader import load_shareseq_multiome
from multiome_agent.eval.grounding import check_grounding
from multiome_agent.eval.shareseq_fault_injection import swap_cell_line_labels
from multiome_agent.eval.shareseq_scoring import score_cell_line_swap

NOVELTY_QUESTION = (
    "This is a pooled multi-cell-line single-cell multiome dataset. Using the QC summary "
    "tool, propose ONE specific, checkable candidate observation or pattern about this "
    "dataset that you find notable -- something with an actual number attached, not a vague "
    "claim. State your confidence."
)


def split_cells_in_half(mdata: MuData, seed: int = 0) -> tuple[list[str], list[str]]:
    """Random 50/50 split of cell barcodes (not stratified -- a plain random
    split is the simplest, least assumption-laden version of CLAUDE.md's
    "split-half replication across random cell halves")."""
    rng = np.random.default_rng(seed)
    idx = rng.permutation(mdata.n_obs)
    names = list(mdata.obs_names)
    half = mdata.n_obs // 2
    return [names[i] for i in idx[:half]], [names[i] for i in idx[half:]]


def _run_novelty_investigation(label: str, mdata_for_tools: MuData, qc_text: str) -> dict:
    result = run_agent(
        NOVELTY_QUESTION, mdata=mdata_for_tools, qc_summary=qc_text, dataset_context=SHARESEQ_DATASET_CONTEXT
    )
    grounding = check_grounding(result)
    return {
        "label": label, "answer": result.answer, "cost_usd": result.estimated_cost_usd,
        "turns": result.turn_count, "grounded": grounding["all_grounded"], "grounding_detail": grounding,
    }


def run_novelty_checks(clean_mdata_with_clusters: MuData, negative_control_swap_fraction: float = 0.6, seed: int = 0) -> dict:
    """Runs the baseline + split-half + negative-control novelty
    investigations for real and returns all results together for direct
    comparison. `clean_mdata_with_clusters` must already have
    `run_shareseq_fixed_core` applied.
    """
    baseline_summary = shareseq_fixed_core_summary(clean_mdata_with_clusters)
    baseline_text = format_shareseq_qc_summary(baseline_summary)

    half1_names, half2_names = split_cells_in_half(clean_mdata_with_clusters, seed=seed)
    half1 = run_shareseq_fixed_core(load_shareseq_multiome()[half1_names].copy())
    half2 = run_shareseq_fixed_core(load_shareseq_multiome()[half2_names].copy())
    half1_text = format_shareseq_qc_summary(shareseq_fixed_core_summary(half1))
    half2_text = format_shareseq_qc_summary(shareseq_fixed_core_summary(half2))

    swapped_mdata, swap_meta = swap_cell_line_labels(
        clean_mdata_with_clusters, fraction=negative_control_swap_fraction, seed=seed
    )
    swap_score = score_cell_line_swap(clean_mdata_with_clusters, swapped_mdata)
    neg_summary = dict(baseline_summary)
    neg_summary["cell_line_recovery_ari_rna"] = swap_score["cell_line_recovery_ari_faulted"]
    neg_text = format_shareseq_qc_summary(neg_summary)

    runs = [
        _run_novelty_investigation("baseline", clean_mdata_with_clusters, baseline_text),
        _run_novelty_investigation("split_half_1", half1, half1_text),
        _run_novelty_investigation("split_half_2", half2, half2_text),
        _run_novelty_investigation("negative_control_cell_line_swap", clean_mdata_with_clusters, neg_text),
    ]
    return {
        "runs": runs,
        "total_cost_usd": sum(r["cost_usd"] for r in runs),
        "negative_control_severity": negative_control_swap_fraction,
        "negative_control_clean_ari": swap_score["cell_line_recovery_ari_clean"],
        "negative_control_faulted_ari": swap_score["cell_line_recovery_ari_faulted"],
    }
