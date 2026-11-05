"""Drug/condition treated-vs-control comparison: a real, tested branch for
a future dataset that has a drug/treatment condition axis (per
`core.condition_detection`) -- neither current dataset (tenx-cell-ranger,
shareseq-multi-cell-lines) actually has one, both are control-only,
confirmed by `detect_condition_groups`, not assumed.

Deliberately reuses `agent.novelty`'s propose+judge machinery end to end
(via `propose_novel_findings`'s `custom_question` parameter and
`judge_novel_findings` unchanged) rather than a parallel pipeline: a
drug-effect finding is still "a specific, checkable pattern in this data,"
just investigated with different tools (`drug_differential_expression`/
`drug_differential_accessibility`/`check_batch_confound`) and a different
question. This also means it's adversarially judged by the SAME skeptical
Judger, whose `_judger_question` was extended to require an independent
`check_batch_confound` call for this specific finding type (see
novelty.py).
"""

from __future__ import annotations

from multiome_agent.agent.loop import AgentRunResult
from multiome_agent.agent.novelty import propose_novel_findings

MAX_CANDIDATES_PER_GROUP_PAIR = 4


def _drug_comparison_question(condition_check: dict) -> str:
    columns = list(condition_check.get("condition_columns_found", {}).keys())
    columns_desc = "\n".join(
        f'- "{col}": real group sizes {sizes}' for col, sizes in condition_check.get("condition_columns_found", {}).items()
    ) or "(none found)"

    return f"""\
This dataset has a real drug/treatment/condition axis -- `check_for_condition_groups` already \
found:

{columns_desc}

For each real condition column above, pick the two most meaningful groups to compare (e.g. a real \
drug/compound arm vs. its real vehicle/DMSO/untreated control -- not two arbitrary doses unless \
that's the only real contrast available). First call `condition_group_qc` on that column to check \
both arms have enough cells; don't build a finding on a near-empty arm.

If this dataset pools multiple cell lines (check `check_for_identity_columns`), you MUST run the \
comparison SEPARATELY within each of a handful of real lines (not pooled across all of them) -- \
ALWAYS pass `cell_line` to `drug_differential_expression`/`drug_differential_accessibility`. A \
pooled treated-vs-control comparison across multiple lines can be ENTIRELY a between-line \
composition confound: if the treated and control arms happen to contain different proportions of \
each real cell line, a gene that's simply a marker of whichever line is over-represented in one arm \
will look like a significant "treatment effect" pooled, with zero real within-line signal at all.

For each real condition column (and, if applicable, each cell line you check it in): call \
`drug_differential_expression` AND `drug_differential_accessibility` for that scope, then call \
`check_batch_confound` on the SAME column/cell_line scope (pass a real batch column from this \
dataset's own metadata if one exists) -- if `likely_confound=true`, SELF-REJECT this candidate \
silently, the differential result could be explained by batch alone. For each candidate gene/peak \
that survives the batch-confound check, run a literature check: `search_pubmed` for the specific \
drug/compound name together with that specific gene (e.g. "{{drug}} {{gene}} expression"), then \
`fetch_pubmed_abstracts` on any promising hit, to check whether the OBSERVED DIRECTION (up or down \
in the treated arm) matches a published drug effect for that gene -- not just that the gene is \
mentioned somewhere near the drug. If the direction DISAGREES with what's published, that's still \
worth recording (as a real, specific disagreement, with `confidence` reflecting how surprising it \
is), not silently discarded -- but if your search finds this EXACT direction already well-reported \
for this drug, SELF-REJECT it as not novel and try a different candidate instead.

Try up to {MAX_CANDIDATES_PER_GROUP_PAIR} distinct candidate genes/peaks this way. Record via \
`record_novel_finding` with `finding_type="positive_relationship"` for a genuine drug-effect claim \
(direction confirmed OR a real, specific disagreement with the literature), with `cell_lines` set \
to every line it was confirmed in, and `evidence` naming the real numbers (logfoldchange, \
pvals_adj, n_treated/n_control, and the batch-confound check's own pvalue/likely_confound). Stop \
once you've recorded 2 findings that survived this way, or once you've tried \
{MAX_CANDIDATES_PER_GROUP_PAIR} distinct candidates, whichever comes first -- ending with fewer \
than 2 (even zero) is a legitimate, expected outcome if every candidate turned out confounded or \
already well-known, not a failure to force past.
"""


def propose_drug_comparison_findings(
    mdata, condition_check: dict, qc_summary: str, dataset_context: str | None = None,
    model: str | None = None, max_turns: int = 30,
) -> tuple[list[dict], AgentRunResult]:
    """The proposal half of the drug-comparison branch -- judging reuses
    `agent.novelty.judge_novel_findings` unchanged (same Judger, extended
    to require `check_batch_confound` for this finding type). Returns
    `([], a zero-cost AgentRunResult)` immediately, with NO agent call at
    all, if `condition_check["is_control_only"]` -- the common case for
    both current datasets, so this branch has zero cost and zero effect
    when there's nothing to investigate.
    """
    if condition_check.get("is_control_only", True):
        return [], AgentRunResult(
            question="(control-only dataset -- no drug/condition axis to investigate)", answer="",
            turn_count=0, hit_max_turns=False, input_tokens=0, output_tokens=0,
            cache_read_tokens=0, cache_creation_tokens=0, estimated_cost_usd=0.0, reasoning_trail_path="",
        )

    return propose_novel_findings(
        mdata, qc_summary, dataset_context=dataset_context, model=model, max_turns=max_turns,
        custom_question=_drug_comparison_question(condition_check),
    )
