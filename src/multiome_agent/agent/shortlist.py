"""Cost-saving shortlist-first path for known-biology checklist generation
and novel-finding pre-filtering.

Motivation: the full `run_agent`-based exploration in `checklist_generator.py`
and `novelty.py` earns its cost by discovering a dataset's identity AND
brainstorming candidate genes/motifs/TF-target pairs for it with no hint --
real work when the identity is genuinely unknown (e.g. the public PBMC
dataset, which carries no cell-line metadata at all). But when a dataset's
real identity is ALREADY retrievable for free from its own metadata
(`core.identity_metadata.detect_identity_columns`, no LLM), most of that
exploration is redundant: a cheap, one-shot recall of well-known markers/TFs
for that identity, verified against THIS dataset's real data via the SAME
already-free deterministic analysis-menu tools, with ONE targeted
literature check per surviving candidate (still free NCBI calls -- just not
an open-ended multi-query brainstorm), reaches the same "confirmed by data"
result far more cheaply.

Every LLM call in this module is a short, tool-minimal, ONE-SHOT call on the
cheapest model (`SHORTLIST_MODEL`, Haiku) -- recall and single-abstract
judgment, never open-ended multi-turn exploration. That is the entire cost
difference versus `agent.loop.run_agent`. The full expensive exploration
remains the fallback whenever this shortlist path doesn't reach enough
confirmed items for a category, or when no free identity is available at
all -- callers never lose coverage, they just pay less for the common case.

Deliberately NOT integrated here: provenance tagging (which items came from
the cheap path vs. the expensive one). The owner explicitly asked for the
report to stay silent on that distinction.
"""

from __future__ import annotations

import anthropic

from multiome_agent.config import get_anthropic_api_key
from multiome_agent.core.clustering import marker_is_recovered
from multiome_agent.core.identity_metadata import detect_identity_columns
from multiome_agent.core.motif_deviations import best_motif_match
from multiome_agent.core.peak_to_gene_links import peak_to_gene_links as _peak_to_gene_links
from multiome_agent.core.regulon_inference import infer_regulon_targets
from multiome_agent.logging_utils import get_logger
from multiome_agent.tools.pubmed import fetch_pubmed_abstracts, search_pubmed

logger = get_logger(__name__)

SHORTLIST_MODEL = "claude-haiku-4-5"

_CATEGORY_INSTRUCTIONS = {
    "rna_marker": "Name RNA marker genes expected to be highly expressed/enriched specifically for this identity.",
    "motif": "Name DNA-binding motifs (identified by the transcription factor that binds them) expected to show elevated chromatin accessibility specifically for this identity.",
    "peak_to_gene": (
        "Name genes reported to be regulated by a specific DISTAL enhancer/regulatory element (NOT "
        "the gene's own promoter), specifically for this identity. `primary` MUST be ONLY the target "
        "gene's bare symbol -- e.g. \"MYC\" -- never a compound description like "
        "\"CELL_LINE_DISTAL_ENHANCER->MYC\" or any string containing \"->\"; the enhancer/distal-element "
        "detail belongs in `claim`, not in `primary`."
    ),
    "regulon_target": (
        "Name specific TF->TARGET gene pairs (a transcription factor reported to directly regulate "
        "ONE specific OTHER named gene's expression -- never the TF's own motif) specifically for "
        "this identity. `primary` MUST be in the EXACT literal format \"TF_SYMBOL->TARGET_SYMBOL\" "
        "with a real \"->\" separator and nothing else -- e.g. \"SPI1->CD14\" or \"GATA3->BCL11B\" -- "
        "never just the TF's bare name alone."
    ),
}

RECORD_SHORTLIST_TOOL = {
    "name": "record_shortlist_candidates",
    "description": (
        "Record your shortlist of well-known candidates from your own background biology "
        "knowledge -- NOT yet verified against any real data, and NOT yet checked against any "
        "literature search you've actually run. Each candidate needs a specific, checkable claim."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "candidates": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "primary": {
                            "type": "string",
                            "description": "Gene symbol, motif/TF name, or 'TF->TARGET' pair (category regulon_target only, e.g. 'SPI1->CD14')",
                        },
                        "claim": {"type": "string", "description": "The specific expected finding, e.g. 'FOXA1 is required for ESR1 chromatin binding in luminal breast cancer'"},
                    },
                    "required": ["primary", "claim"],
                },
            },
        },
        "required": ["candidates"],
    },
}

JUDGE_ABSTRACT_TOOL = {
    "name": "record_abstract_judgment",
    "description": "Record whether the abstract actually, specifically supports the claim.",
    "input_schema": {
        "type": "object",
        "properties": {"supports_claim": {"type": "boolean"}},
        "required": ["supports_claim"],
    },
}


def resolve_cheap_identity(mdata) -> list[str] | None:
    """Real cell-line/identity names already resolvable for free from this
    dataset's own metadata (`detect_identity_columns`, no LLM) -- returns
    distinct identities ordered by cell count (most abundant first), or
    None if nothing usable was found (e.g. the public PBMC dataset carries
    no cell-line metadata at all -- callers should fall back to the full
    agent-exploration path unconditionally in that case) or if `mdata`
    doesn't have the expected structure (defensive: this is a cost
    optimization, it must never be the reason the pipeline crashes).
    """
    try:
        found = detect_identity_columns(mdata).get("identity_columns_found", {})
    except Exception:
        logger.debug("resolve_cheap_identity: detect_identity_columns failed, treating as unavailable", exc_info=True)
        return None

    for col_name, value_counts in found.items():
        if col_name.endswith("cell_line_name"):
            names = {v: c for v, c in value_counts.items() if not str(v).upper().startswith("ACH-")}
            if names:
                return sorted(names, key=lambda n: -names[n])
    for col_name, value_counts in found.items():
        if col_name.endswith("may_have_wrong_cell_line_label_based_on_rna_cluster"):
            continue
        sample = next(iter(value_counts), None)
        if sample is not None and not str(sample).replace(".", "").replace("-", "").isdigit():
            return sorted(value_counts, key=lambda n: -value_counts[n])
    return None


def _call_shortlist_model(prompt: str, tool: dict, model: str) -> dict | None:
    client = anthropic.Anthropic(api_key=get_anthropic_api_key())
    response = client.messages.create(
        model=model, max_tokens=1024, tools=[tool],
        tool_choice={"type": "tool", "name": tool["name"]},
        messages=[{"role": "user", "content": prompt}],
    )
    for block in response.content:
        if block.type == "tool_use" and block.name == tool["name"]:
            return block.input
    return None


def generate_shortlist_candidates(identity: str, category: str, n: int = 5, model: str = SHORTLIST_MODEL) -> list[dict]:
    """One-shot, tool-minimal recall (no exploration, no dataset access) of
    up to `n` well-known candidates for `category` in the context of
    `identity`. NOT yet grounded in this dataset or in a real literature
    search -- `shortlist_checklist_items_for_category` below verifies each
    one against both before it's ever recorded.
    """
    prompt = (
        f"From your own background biology knowledge (do not search anything, just recall): "
        f"for the cell line/identity \"{identity}\", {_CATEGORY_INSTRUCTIONS[category]} "
        f"List up to {n} distinct, specific, well-established candidates, most confident first."
    )
    result = _call_shortlist_model(prompt, RECORD_SHORTLIST_TOOL, model)
    candidates = (result or {}).get("candidates", [])
    logger.info("generate_shortlist_candidates: %s/%s -> %d candidates", identity, category, len(candidates))
    return candidates[:n]


def judge_abstract_supports_claim(claim: str, abstract: str, model: str = SHORTLIST_MODEL) -> bool:
    """One-shot, single-abstract judgment -- NOT a substitute for the full
    adversarial Judger; just answers "does this real abstract's text
    actually support this specific claim," the same question the expensive
    exploration asks per candidate, bounded here to one abstract instead of
    an open-ended multi-query search.
    """
    prompt = (
        f"Claim: {claim}\n\nAbstract:\n{abstract}\n\n"
        "Does this abstract actually, specifically support the claim above (not just a loosely "
        "related topic, and not just because it mentions the same gene)? Record your judgment."
    )
    result = _call_shortlist_model(prompt, JUDGE_ABSTRACT_TOOL, model)
    return bool((result or {}).get("supports_claim"))


def verify_candidate_in_data(mdata, category: str, primary: str) -> bool:
    """Deterministic, zero-LLM check of whether `primary` actually shows the
    expected signal in THIS dataset -- reuses the same already-free
    analysis-menu tools the expensive agent path would call, just invoked
    directly instead of through a model's tool-use turn.
    """
    try:
        if category == "rna_marker":
            return marker_is_recovered(mdata.mod["rna"], primary.strip().upper())
        if category == "motif":
            names = mdata.mod["atac"].uns.get("chromvar_motif_names")
            return bool(names) and best_motif_match(list(names), primary) is not None
        if category == "peak_to_gene":
            result = _peak_to_gene_links(mdata, primary)
            return "error" not in result and result["any_significant_distal_link"]
        if category == "regulon_target":
            tf, _, target = primary.partition("->")
            if not target:
                return False
            target = target.strip().upper()
            # `target_gene=target` guarantees this specific gene's real entry
            # appears in the (otherwise capped) `targets` list -- see
            # `infer_regulon_targets`'s docstring for why that cap exists.
            result = infer_regulon_targets(mdata, tf.strip(), target_gene=target)
            if "error" in result:
                return False
            return any(t["gene"] == target and t["significant"] for t in result["targets"])
    except Exception:
        logger.exception("verify_candidate_in_data failed for %s/%s", category, primary)
    return False


def shortlist_checklist_items_for_category(mdata, identity: str, category: str, needed: int = 3, max_candidates: int = 5, model: str = SHORTLIST_MODEL):
    """The full cheap path for one checklist category: recall -> verify in
    data (free) -> ONE targeted literature check per data-verified candidate
    (free search + one cheap judgment) -> up to `needed` confirmed
    `ChecklistItem`s. Returns fewer than `needed` (even zero) if the
    shortlist doesn't have enough real, confirmable candidates -- that's an
    honest result the caller falls back to expensive exploration to fill,
    not a bug.
    """
    from multiome_agent.agent.checklist_generator import ChecklistItem  # local import avoids a module cycle

    confirmed = []
    for candidate in generate_shortlist_candidates(identity, category, n=max_candidates, model=model):
        if len(confirmed) >= needed:
            break
        primary, claim = candidate.get("primary", ""), candidate.get("claim", "")
        if not primary or not claim:
            continue
        # Defensive normalization against a real observed failure mode: despite the
        # prompt's explicit format requirement, a cheap model can still return a
        # compound string for "peak_to_gene" (e.g. "CELL_LINE_DISTAL_ENHANCER->MYC"
        # instead of bare "MYC") -- recover the real gene symbol (text after the
        # LAST "->") rather than let a free-but-wasted shortlist attempt silently
        # fall through to the expensive fallback path for a trivially fixable format
        # slip. "regulon_target" has no equivalent safe recovery (there's no way to
        # guess a missing target from a bare TF name), so a missing "->" there is
        # left to fail honestly via `verify_candidate_in_data`.
        if category == "peak_to_gene" and "->" in primary:
            primary = primary.rsplit("->", 1)[-1].strip()
        if not verify_candidate_in_data(mdata, category, primary):
            continue
        hits = search_pubmed(claim)
        if not hits:
            continue
        pmid = hits[0]["pmid"]
        abstract_data = fetch_pubmed_abstracts([pmid]).get(pmid)
        if not abstract_data or not abstract_data.get("abstract"):
            continue
        if not judge_abstract_supports_claim(claim, abstract_data["abstract"], model=model):
            continue
        confirmed.append(ChecklistItem(
            category=category, gene_or_motif=primary, claim=claim, pmid=pmid,
            journal=abstract_data.get("journal", ""), year=str(abstract_data.get("year", "")),
            confirmed_present_in_data=True,
        ))
    logger.info(
        "shortlist_checklist_items_for_category: %s/%s -> %d/%d confirmed from shortlist",
        identity, category, len(confirmed), needed,
    )
    return confirmed


def cheap_already_known_check(finding: str, model: str = SHORTLIST_MODEL) -> dict | None:
    """ONE free `search_pubmed` call + (if it hits) one cheap single-abstract
    judgment -- catches a proposed "novel" finding that's actually trivially
    already well-established, without spending the full adversarial Judger
    agent call just to rediscover that. Returns a struck-down verdict dict
    (same shape `judge_novel_findings` expects) if the cheap check finds
    it's already known; None if unresolved, meaning the caller should
    escalate to the real Judger.
    """
    hits = search_pubmed(finding)
    if not hits:
        return None
    pmid = hits[0]["pmid"]
    abstract_data = fetch_pubmed_abstracts([pmid]).get(pmid)
    if not abstract_data or not abstract_data.get("abstract"):
        return None
    if not judge_abstract_supports_claim(finding, abstract_data["abstract"], model=model):
        return None
    return {
        "likely_artifact": False,
        "artifact_reasoning": "(not evaluated -- struck down by the cheap already-known prefilter before any artifact-checking)",
        "already_known": True,
        "novelty_reasoning": (
            f"Cheap prefilter found PMID {pmid} ({abstract_data.get('journal', '')}, "
            f"{abstract_data.get('year', '')}) already describing this."
        ),
        "verdict": "struck_down",
    }
