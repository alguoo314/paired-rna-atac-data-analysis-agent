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
from multiome_agent.agent.shortlist import (
    SHORTLIST_MODEL,
    resolve_cheap_identity,
    shortlist_checklist_items_for_category,
)
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

CHECKLIST_CATEGORIES = ("rna_marker", "motif", "peak_to_gene", "regulon_target")

RECORD_CHECKLIST_ITEM_TOOL = {
    "name": "record_checklist_item",
    "description": (
        "Record ONE known-biology checklist item, grounded in a literature abstract you "
        "actually fetched and read (via fetch_pubmed_abstracts). Two cases: (1) a claim you "
        "verified holds in this dataset -- confirmed_present_in_data=true; (2) a specific, "
        "real candidate whose own data-verification tool call CLEARLY contradicted the "
        "literature claim (not significant, wrong sign, or the gene/motif absent from this "
        "dataset) -- confirmed_present_in_data=false. Do NOT record a candidate that checked "
        "out fine but simply wasn't needed once enough confirmed items were found in that "
        "category -- that's not a rejection, just leave it unrecorded. Call this once per item, "
        "up to 4 times per category (matching the 4-candidate-per-category cap)."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "category": {"type": "string", "enum": list(CHECKLIST_CATEGORIES)},
            "gene_or_motif": {
                "type": "string",
                "description": (
                    "Gene symbol or motif name this item is about. For category \"regulon_target\" "
                    "specifically, use \"TF->TARGET\" (e.g. \"SPI1->CD14\") since that item is about "
                    "a specific TF-target gene pair, not one gene alone."
                ),
            },
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
distinct cell lines/lineages, don't try to exhaustively characterize every one, and don't \
concentrate the whole checklist on a single line either -- identify the overall composition \
quickly, then span your candidates across a handful of different lines chosen at random rather \
than defaulting to whichever one is easiest or most already-characterized. Spanning means the SET \
of claims should name different lines across it -- e.g. one claim about line A, the next about \
line B, the next about line C -- NOT a single claim that vaguely covers several lines at once \
(e.g. "elevated in epithelial cancer cells" is too generic; "ESR1 is expressed in T-47D, an \
ERalpha-positive luminal breast cancer line" is the right level of specificity). Every individual \
claim should be exactly as specific and well-established as if you'd focused on one line the \
whole time -- spanning changes WHICH line each claim is about, not how specific any one claim is.

Then, for that specific identity, build a known-biology checklist grounded in real literature, \
working through these four categories IN ORDER:

1. "rna_marker" -- an RNA marker gene expected to be expressed/enriched for this specific identity.
2. "motif" -- a DNA-binding motif expected to show elevated chromatin accessibility for this \
specific identity.
3. "peak_to_gene" -- a specific, literature-reported distal regulatory element (an enhancer or \
other distal regulatory region, not the gene's own promoter) reported to regulate a named gene's \
expression for this specific identity. Verify with `peak_to_gene_links`: a real link requires \
`any_significant_distal_link=true` for that gene, i.e. at least one distal peak's accessibility \
actually tracks that gene's own expression in THIS dataset -- not just that the gene itself is \
expressed.
4. "regulon_target" -- a specific, literature-reported TF-target-gene pair (e.g. "SPI1 regulates \
CD14 in myeloid cells") for this specific identity -- NOT the TF's own motif, a DIFFERENT gene it's \
reported to regulate. Verify with `regulon_inference(tf_gene=TF, target_gene=TARGET)` -- ALWAYS \
pass `target_gene` once you have a specific candidate in mind, not just `tf_gene` alone: the \
`targets` list is capped (an abundant TF motif can have thousands of real candidates), and \
`target_gene` guarantees your specific candidate's real entry appears in the result even if it \
wouldn't otherwise make that cap. Require `significant=true` for it in the result -- i.e. the \
TF's own RNA actually tracks that SPECIFIC target gene's own RNA in THIS dataset, with real \
motif/chromatin evidence supporting candidacy, not just that both genes happen to be expressed.

Deliberately NOT a category here: a transcription factor's own RNA tracking its own motif's \
accessibility (`tf_motif_correlation`). That used to be category 3, but it's a narrower, \
easily-already-known claim compared to "regulon_target" above (a TF tracking a SPECIFIC OTHER \
gene it regulates) -- `tf_motif_correlation` remains a real, usable tool for other purposes (e.g. \
corroborating a candidate), just not something to build a checklist category around anymore.

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
(`peak_to_gene_links`, `regulon_inference`, \
`cross_modal_marker_check`, or `top_cluster_markers`) before recording; if the gene/motif turns \
out absent, pick a different candidate rather than forcing it. Don't spend excessive turns \
re-verifying a candidate that already checked out, and don't backtrack to a finished category \
once you've moved on -- but do move through all 4 categories with real effort on each, not just \
the first one you reach.

Report failures, not just successes: if a specific candidate's own data-verification tool call \
clearly CONTRADICTS the literature claim -- not significant, wrong sign, or absent from this \
dataset -- record that too via `record_checklist_item` with `confirmed_present_in_data=false`, \
so a real negative result is captured rather than silently discarded. This is different from a \
candidate that simply checked out fine but wasn't needed once you'd already found enough \
confirmed items in that category -- don't record those, just move on.
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


_FALLBACK_CATEGORY_DESCRIPTIONS = {
    "rna_marker": "an RNA marker gene expected to be expressed/enriched for this specific identity.",
    "motif": "a DNA-binding motif expected to show elevated chromatin accessibility for this specific identity.",
    "peak_to_gene": (
        "a specific, literature-reported distal regulatory element (not the gene's own promoter) "
        "reported to regulate a named gene's expression for this specific identity. Verify with "
        "`peak_to_gene_links`: requires `any_significant_distal_link=true` for that gene."
    ),
    "regulon_target": (
        "a specific, literature-reported TF-target-gene pair (NOT the TF's own motif) for this "
        "specific identity. Verify with `regulon_inference(tf_gene=TF, target_gene=TARGET)` -- "
        "ALWAYS pass `target_gene` so your specific candidate is guaranteed to appear in the "
        "(capped) `targets` list -- and require `significant=true` for it."
    ),
}


def _fallback_question(categories_needed: dict[str, int], identity_hint: str | None) -> str:
    """Builds a CHECKLIST_QUESTION-equivalent scoped to ONLY the categories
    that still need filling after the cheap shortlist pass -- see
    `generate_checklist`'s docstring. Passing `identity_hint` (already
    resolved for free from this dataset's own metadata) skips re-deriving
    identity from scratch, saving the turns that would otherwise go into
    `check_for_identity_columns`/marker-based inference.
    """
    lines = [f'- "{cat}": {_FALLBACK_CATEGORY_DESCRIPTIONS[cat]} Still need {n} more well-grounded item(s) in this category.' for cat, n in categories_needed.items()]
    identity_line = (
        f'This dataset\'s real identity has already been determined from its own metadata: '
        f'"{identity_hint}". Use it directly -- no need to re-derive it.'
        if identity_hint else
        "First, determine this dataset's real identity from evidence (check_for_identity_columns "
        "FIRST, then marker-based inference if that finds nothing) -- never told directly."
    )
    return f"""\
{identity_line}

Build a known-biology checklist grounded in real literature for ONLY the following categories \
(other categories were already covered by a cheaper pass and need no further attention from you):

{chr(10).join(lines)}

For each candidate: `search_pubmed`, then `fetch_pubmed_abstracts` on the promising PMIDs -- a \
title alone is never sufficient grounding, only record a claim backed by a real, read abstract. \
Verify data-presence with a real tool call (`tf_motif_correlation`, `peak_to_gene_links`, \
`regulon_inference`, `cross_modal_marker_check`, or `top_cluster_markers`) before recording; if \
the gene/motif turns out absent, pick a different candidate rather than forcing it. Try up to 4 \
distinct candidates per category before giving up on filling it -- that's a real, honest result, \
not a failure. Report failures too: if a specific candidate's own data-verification tool call \
clearly CONTRADICTS the literature claim, record that via `record_checklist_item` with \
`confirmed_present_in_data=false`.
"""


def generate_checklist(
    mdata, qc_summary: str, dataset_context: str | None = None, model: str | None = None, max_turns: int = 60,
    categories: tuple[str, ...] = CHECKLIST_CATEGORIES, existing_items: list[ChecklistItem] | None = None,
    use_shortlist: bool = True, shortlist_model: str = SHORTLIST_MODEL,
) -> tuple[list[ChecklistItem], AgentRunResult]:
    """`dataset_context` defaults to `OWN_DATA_CONTEXT` -- the right choice for
    arbitrary data, since it lets the agent name a real identity it determines
    (e.g. a cell line) while still forbidding it from ever naming the data's
    source. Callers running the PBMC report pass `PBMC_DATASET_CONTEXT`
    explicitly instead, since that dataset has no cell-line identity to name.

    `categories` -- which categories to (re)generate; defaults to all five.
    Pass a subset (e.g. `("peak_to_gene", "regulon_target")`) together with
    `existing_items` (the previously confirmed items for the OTHER
    categories, e.g. loaded from a prior report run) to add new categories
    to an existing checklist WITHOUT re-running or re-paying for categories
    that haven't changed -- those are returned completely untouched,
    verbatim, never re-verified or re-cited.

    `use_shortlist` -- if True (default) and this dataset's identity is
    cheaply resolvable from its own metadata (`agent.shortlist.
    resolve_cheap_identity`, free, no LLM), tries a cheap shortlist-first
    path per category (one-shot recall on `shortlist_model` + free data
    verification via the existing deterministic analysis-menu tools + one
    targeted literature check) before falling back to the full open-ended
    `run_agent` exploration for whatever a category's shortlist pass didn't
    fill. If no cheap identity is available at all (e.g. the public PBMC
    dataset), this has no effect -- behavior is identical to before this
    existed. Either way, every category still ends up with as many
    confirmed items as the expensive path alone would have found; the
    shortlist only ever reduces cost, never coverage.
    """
    preserved = [i for i in (existing_items or []) if i.category not in categories]

    shortlist_items: list[ChecklistItem] = []
    remaining_needed: dict[str, int] = {c: 3 for c in categories}
    identity_hint = None
    if use_shortlist:
        identities = resolve_cheap_identity(mdata)
        if identities:
            identity_hint = identities[0]
            for category in categories:
                found = shortlist_checklist_items_for_category(mdata, identity_hint, category, needed=3, model=shortlist_model)
                shortlist_items.extend(found)
                remaining_needed[category] = max(0, 3 - len(found))

    categories_needing_fallback = {c: n for c, n in remaining_needed.items() if n > 0}
    if categories_needing_fallback:
        question = _fallback_question(categories_needing_fallback, identity_hint)
        result = run_agent(
            question, model=model, max_turns=max_turns, mdata=mdata, qc_summary=qc_summary,
            tools=CHECKLIST_TOOLS, dataset_context=dataset_context if dataset_context is not None else OWN_DATA_CONTEXT,
        )
        fallback_items = [
            ChecklistItem(**tc["input"])
            for tc in result.tool_calls
            if tc["name"] == "record_checklist_item" and not tc.get("is_error", False)
        ]
    else:
        fallback_items = []
        result = AgentRunResult(
            question="(shortlist-only -- no agent exploration needed for any category)", answer="",
            turn_count=0, hit_max_turns=False, input_tokens=0, output_tokens=0,
            cache_read_tokens=0, cache_creation_tokens=0, estimated_cost_usd=0.0, reasoning_trail_path="",
        )

    # Deterministic dedup + per-category cap, not relied on either source's
    # own discipline: a real run once recorded two "tf_motif_tracking" items
    # both about FOXA1 (two different literature framings of the same gene)
    # -- keep only the first item per (category, gene) pair. Cross-category
    # reuse of the same gene is fine and expected (e.g. FOXA1 as both a
    # "motif" item and a "tf_motif_tracking" item are two distinct,
    # non-redundant claims), so dedup is scoped to (category, gene), not
    # gene alone. The cap matters now specifically because two independent
    # sources (shortlist + fallback) can both contribute to the same
    # category.
    seen: set[tuple[str, str]] = set()
    per_category_count: dict[str, int] = {}
    items = []
    for item in shortlist_items + fallback_items:
        key = (item.category, item.gene_or_motif.lower())
        if key in seen or per_category_count.get(item.category, 0) >= 3:
            continue
        seen.add(key)
        per_category_count[item.category] = per_category_count.get(item.category, 0) + 1
        items.append(item)

    logger.info(
        "Checklist generated: %d new items (%s) [+%d preserved from existing_items], cost=$%.5f",
        len(items), ", ".join(f"{c}={sum(1 for i in items if i.category == c)}" for c in categories),
        len(preserved), result.estimated_cost_usd,
    )
    return preserved + items, result
