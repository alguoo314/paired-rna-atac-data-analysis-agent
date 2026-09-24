"""Known-biology checklist generation via literature RAG.

Per the resume-project rebuild: the agent determines this dataset's real
identity itself (see agent.prompts principle 6 / core.cross_modal_validation
-- it is never told directly), then finds real, actually-read (not just
title-matched) literature specific to that identity and records up to 3
checklist items per category as STRUCTURED tool calls, not free text, so
extraction afterward is exact rather than a regex over prose (same pattern
as agent.loader_selector's decide_loading_strategy).

This deliberately reuses `agent.loop.run_agent` (its existing tool-execution
machinery, cost/turn tracking, reasoning trail) rather than a bespoke loop --
unlike `loader_selector`, this step genuinely needs the full analysis
toolset (marker discovery, cross-modal validation, motif correlation) plus
literature RAG, which `run_agent` already wires up.
"""

from __future__ import annotations

from dataclasses import dataclass

from multiome_agent.agent.loop import TOOLS, AgentRunResult, run_agent
from multiome_agent.agent.prompts import OWN_DATA_CONTEXT
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

CHECKLIST_CATEGORIES = ("rna_marker", "motif", "tf_motif_tracking")

RECORD_CHECKLIST_ITEM_TOOL = {
    "name": "record_checklist_item",
    "description": (
        "Record ONE known-biology checklist item, grounded in a literature abstract you "
        "actually fetched and read (via fetch_pubmed_abstracts) and verified is usable in "
        "THIS dataset. Call this once per item, up to 3 times per category."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "category": {"type": "string", "enum": list(CHECKLIST_CATEGORIES)},
            "gene_or_motif": {"type": "string", "description": "Gene symbol or motif name this item is about"},
            "claim": {"type": "string", "description": "The specific expected finding, e.g. 'ESR1 is highly expressed in ER+ breast lineage cells'"},
            "pmid": {"type": "string"},
            "journal": {"type": "string"},
            "year": {"type": "string"},
            "confirmed_present_in_data": {
                "type": "boolean",
                "description": "Whether you verified this gene/motif is actually usable in this dataset (via a real tool call), not just assumed",
            },
        },
        "required": ["category", "gene_or_motif", "claim", "pmid", "journal", "year", "confirmed_present_in_data"],
    },
}

CHECKLIST_TOOLS = TOOLS + [RECORD_CHECKLIST_ITEM_TOOL]

CHECKLIST_QUESTION = """\
First, determine this dataset's real identity. Check `check_for_identity_columns` FIRST -- if \
this dataset's own metadata already encodes identity (e.g. a real cell-line column, or a column \
of DepMap "ACH-XXXXXX" IDs you resolve with `resolve_depmap_id`), that's real discovery and you \
should use it directly, not re-derive it the hard way. Only fall back to marker-based inference \
(cluster markers, cross-modal validation) if no such field exists. If this dataset pools multiple \
distinct cell lines/lineages, don't try to exhaustively characterize every one -- identify the \
overall composition quickly, then pick ONE clearly-characterized line to build the checklist \
around.

Then, for that specific identity, build a known-biology checklist grounded in real literature, \
working through these three categories IN ORDER:

1. "rna_marker" -- an RNA marker gene expected to be expressed/enriched for this specific identity.
2. "motif" -- a DNA-binding motif expected to show elevated chromatin accessibility for this \
specific identity.
3. "tf_motif_tracking" -- a transcription factor whose own RNA expression is reported (ideally \
in the same paper) to track its own motif's accessibility, for this specific identity.

AIM FOR 3 well-grounded items in EACH category -- that's the real goal, not a stretch target. \
Only settle for fewer after a genuine effort: try a HANDFUL of different candidate genes/motifs \
per category (different search queries, not just rephrasing the same one) before concluding a \
category can't reach 3. But this effort is BOUNDED, not open-ended: at most 4 distinct candidates \
per category. If none of your first 4 distinct candidates yields a groundable, data-present item, \
STOP and move on with whatever you have -- don't keep expanding the search for a 5th, 6th candidate, \
that's exactly the failure mode that once burned $6 and 60 turns on a single category. Never record \
two items in the same category about the same gene/motif -- if your best candidate already has an \
item recorded, that category is done regardless of how many turns remain; move to the next one. \
It's fine to end up with fewer than 3 in a category once you've hit the 4-candidate cap -- that's \
a real, honest result, not a failure. For each candidate claim: `search_pubmed`, then \
`fetch_pubmed_abstracts` on the promising PMIDs -- a title alone is never sufficient grounding, \
only record a claim backed by a real, read abstract. Verify data-presence with a real tool call \
(`tf_motif_correlation`, \
`cross_modal_marker_check`, or `top_cluster_markers`) before recording; if the gene/motif turns \
out absent, pick a different candidate rather than forcing it. Don't spend excessive turns \
re-verifying a candidate that already checked out, and don't backtrack to a finished category \
once you've moved on -- but do move through all 3 categories with real effort on each, not just \
the first one you reach.
"""


@dataclass
class ChecklistItem:
    category: str
    gene_or_motif: str
    claim: str
    pmid: str
    journal: str
    year: str
    confirmed_present_in_data: bool


def generate_checklist(
    mdata, qc_summary: str, dataset_context: str | None = None, model: str | None = None, max_turns: int = 60
) -> tuple[list[ChecklistItem], AgentRunResult]:
    """`dataset_context` defaults to `OWN_DATA_CONTEXT` -- the right choice for
    arbitrary data, since it lets the agent name a real identity it determines
    (e.g. a cell line) while still forbidding it from ever naming the data's
    source. Callers running the PBMC report pass `PBMC_DATASET_CONTEXT`
    explicitly instead, since that dataset has no cell-line identity to name.
    """
    result = run_agent(
        CHECKLIST_QUESTION, model=model, max_turns=max_turns, mdata=mdata, qc_summary=qc_summary,
        tools=CHECKLIST_TOOLS, dataset_context=dataset_context if dataset_context is not None else OWN_DATA_CONTEXT,
    )
    all_items = [
        ChecklistItem(**tc["input"])
        for tc in result.tool_calls
        if tc["name"] == "record_checklist_item" and not tc.get("is_error", False)
    ]
    # Deterministic dedup, not relied on the model's own discipline: a real
    # run recorded two "tf_motif_tracking" items both about FOXA1 (two
    # different literature framings of the same gene) while still hitting
    # max_turns searching for a genuinely different 3rd candidate -- keep
    # only the first item per (category, gene) pair rather than trust the
    # agent never submits a within-category duplicate. Cross-category reuse
    # of the same gene is fine and expected (e.g. FOXA1 as both a "motif"
    # item and a "tf_motif_tracking" item are two distinct, non-redundant
    # claims), so dedup is scoped to (category, gene), not gene alone.
    seen: set[tuple[str, str]] = set()
    items = []
    for item in all_items:
        key = (item.category, item.gene_or_motif.lower())
        if key in seen:
            continue
        seen.add(key)
        items.append(item)
    logger.info(
        "Checklist generated: %d items (%s), cost=$%.5f",
        len(items), ", ".join(f"{c}={sum(1 for i in items if i.category == c)}" for c in CHECKLIST_CATEGORIES),
        result.estimated_cost_usd,
    )
    return items, result
