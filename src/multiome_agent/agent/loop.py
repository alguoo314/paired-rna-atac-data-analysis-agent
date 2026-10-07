"""Minimal hand-rolled tool-calling agent loop (step 5).

Hand-rolled rather than the SDK's beta `tool_runner` helper: this project
needs a precise custom max-turns guard and a structured, readable
reasoning-trail log (CLAUDE.md's "everything is logged" principle), both of
which a manual `while` loop over `client.messages.create` makes
straightforward to instrument exactly where needed, without taking on a
beta SDK dependency for something this small.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import anthropic

from multiome_agent.agent.fixed_core_cache import get_agent_fixed_core
from multiome_agent.agent.prompts import CORE_SYSTEM_PROMPT, PBMC_DATASET_CONTEXT
from multiome_agent.agent.qc_summary import fixed_core_summary, format_qc_summary
from multiome_agent.config import AGENT_MODEL, MAX_TURNS, REPO_ROOT, get_anthropic_api_key
from multiome_agent.core.clustering import ATAC_CLUSTER_KEY, RNA_CLUSTER_KEY, top_gene_activity_cluster_markers, top_rna_cluster_markers
from multiome_agent.core.condition_detection import condition_group_qc, detect_condition_groups
from multiome_agent.core.cross_modal_validation import cross_modal_marker_validation
from multiome_agent.core.identity_metadata import detect_identity_columns
from multiome_agent.logging_utils import get_logger
from multiome_agent.menu.peak_to_gene_links import peak_to_gene_links
from multiome_agent.menu.regulon_inference import regulon_inference
from multiome_agent.menu.tf_motif_correlation import tf_motif_correlation
from multiome_agent.tools.code_execution import run_analysis_code
from multiome_agent.tools.depmap import resolve_depmap_id
from multiome_agent.tools.enrichr import enrich_gene_set
from multiome_agent.tools.pubmed import fetch_pubmed_abstracts, search_pubmed

logger = get_logger(__name__)

RUNS_DIR = REPO_ROOT / "logs" / "agent_runs"

# $ per 1M tokens, verified current rates (not training-data memory) for
# step 10's model comparison table.
MODEL_COSTS = {
    "claude-haiku-4-5": {"input": 1.00, "output": 5.00},
    "claude-sonnet-5": {"input": 3.00, "output": 15.00},
    "claude-opus-5": {"input": 5.00, "output": 25.00},
    "claude-fable-5": {"input": 10.00, "output": 50.00},
}

TOOLS = [
    {
        "name": "tf_motif_correlation",
        "description": (
            "Correlate a transcription factor's own RNA expression with its own "
            "DNA-binding motif's chromVAR accessibility deviation score, across "
            "cells in this PBMC multiome dataset. Returns spearman_rho and pvalue, "
            "or an error if the gene isn't in the data or has no matching motif."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {"type": "string", "description": "TF gene symbol, e.g. SPI1"},
            },
            "required": ["gene"],
        },
    },
    {
        "name": "peak_to_gene_links",
        "description": (
            "Test whether any DISTAL peak's accessibility (NOT the gene's own gene-activity "
            "score, which only covers the gene body + promoter) tracks a gene's own RNA "
            "expression, across cells in this dataset. Use this to actually check a "
            "distal-enhancer-regulation explanation for RNA-ATAC discordance, rather than just "
            "naming it as a possibility. Returns the top candidate linked peaks (distance, "
            "Spearman rho, BH-corrected q-value) and whether any cleared BOTH a q<0.05 AND an "
            "|rho|>=0.2 bar (q-value alone is too permissive at real cell counts), or an "
            "error if the gene isn't in the data or has no GENCODE coordinates."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {"type": "string", "description": "Gene symbol, e.g. CD14"},
            },
            "required": ["gene"],
        },
    },
    {
        "name": "regulon_inference",
        "description": (
            "Given a transcription factor gene symbol, finds its candidate target genes via real "
            "motif + chromatin evidence (the TF's own JASPAR motif present in a target's own "
            "promoter, or in a distal peak that itself significantly links to that target's "
            "expression), then reports, PER TARGET GENE, the Spearman correlation between the "
            "TF's own RNA and that SPECIFIC target gene's own RNA -- never an aggregate score "
            "across many genes, and never the TF's own expression standing in for a target. Use "
            "this for a structurally different, harder-to-already-know claim than "
            "tf_motif_correlation: does this TF's RNA actually track a SPECIFIC other gene it's "
            "predicted to regulate? The returned `targets` list is CAPPED (an abundant TF motif "
            "can have thousands of real candidates) -- if you already have one specific target "
            "gene in mind (e.g. verifying a literature-reported TF-target pair), pass "
            "`target_gene` so its real entry is guaranteed to appear even if it wouldn't "
            "otherwise make the cap. Returns an error if the TF isn't in the data, has no matching "
            "motif, or this dataset's motif-match data isn't available."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "tf_gene": {"type": "string", "description": "Transcription factor gene symbol, e.g. SPI1"},
                "target_gene": {
                    "type": "string",
                    "description": "Optional: a specific target gene you want a guaranteed answer for, e.g. CD14",
                },
            },
            "required": ["tf_gene"],
        },
    },
    {
        "name": "search_pubmed",
        "description": (
            "Search PubMed for articles relevant to a query. Returns up to 5 hits (title, "
            "journal, year, PMID) -- TITLES ONLY, for finding candidates. A title match alone "
            "is not evidence for a claim; call `fetch_pubmed_abstracts` on the promising PMIDs "
            "before citing one, so the claim is grounded in what the paper actually says."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {"type": "string", "description": "PubMed search query"},
            },
            "required": ["query"],
        },
    },
    {
        "name": "fetch_pubmed_abstracts",
        "description": (
            "Fetch the real abstract text for a short list of PMIDs (from a prior "
            "`search_pubmed` call). Read the actual content before citing a paper as "
            "supporting a specific claim -- a title alone can be misleading about what a "
            "paper actually found."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "pmids": {"type": "array", "items": {"type": "string"}, "description": "PMIDs to fetch abstracts for"},
            },
            "required": ["pmids"],
        },
    },
    {
        "name": "enrich_gene_set",
        "description": (
            "Run Gene Ontology Biological Process enrichment (via Enrichr) on a list "
            "of gene symbols. Returns up to 10 enriched terms with adjusted p-values."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "genes": {"type": "array", "items": {"type": "string"}, "description": "Gene symbols"},
            },
            "required": ["genes"],
        },
    },
    {
        "name": "resolve_depmap_id",
        "description": (
            "Resolve a Broad Institute DepMap Model ID (format \"ACH-XXXXXX\") to its real cell-line "
            "name and disease/lineage annotation, via a live lookup against EBI's public Cellosaurus "
            "cell-line database. Use this if you notice a metadata column containing ACH-formatted "
            "values."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"depmap_id": {"type": "string", "description": "e.g. 'ACH-000001'"}},
            "required": ["depmap_id"],
        },
    },
    {
        "name": "check_for_identity_columns",
        "description": (
            "Check whether this dataset's own metadata already carries an identity-encoding column "
            "(e.g. a real cell-line/lineage/genotype/donor field) -- by column name, then real values "
            "if found. Call this FIRST, before doing exhaustive marker-based cluster analysis: if this "
            "dataset already tells you its own identity through a real data field, using that is "
            "legitimate discovery (you found it in the data, you weren't told it), not a shortcut "
            "around figuring it out. If nothing is found, fall back to marker-based inference "
            "(top_cluster_markers, cross_modal_marker_check)."
        ),
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "list_clusters",
        "description": (
            "List the real RNA and ATAC Leiden cluster IDs for this dataset -- call this "
            "before top_cluster_markers/top_gene_activity_markers/cross_modal_marker_check "
            "so you pass a cluster ID that actually exists."
        ),
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "top_cluster_markers",
        "description": (
            "Top RNA marker genes of one Leiden RNA cluster, by differential-expression "
            "score (gene, log-fold-change, adjusted p-value). Use this to figure out what a "
            "cluster actually IS from the data itself, rather than assuming -- this is how "
            "you determine cell type / cell line / condition identity when it hasn't been "
            "told to you directly."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "cluster": {"type": "string", "description": "RNA Leiden cluster ID, from list_clusters"},
                "n": {"type": "integer", "description": "How many top markers to return (default 15)"},
            },
            "required": ["cluster"],
        },
    },
    {
        "name": "top_gene_activity_markers",
        "description": (
            "Top ATAC gene-activity markers of one Leiden ATAC cluster (same idea as "
            "top_cluster_markers, computed on chromatin accessibility near each gene instead "
            "of RNA expression) -- the independent, ATAC-side signal for the same identity "
            "question."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "cluster": {"type": "string", "description": "ATAC Leiden cluster ID, from list_clusters"},
                "n": {"type": "integer", "description": "How many top markers to return (default 15)"},
            },
            "required": ["cluster"],
        },
    },
    {
        "name": "cross_modal_marker_check",
        "description": (
            "Textbook ArchR/Signac cross-modal cell-type-call validation for one gene: is it "
            "a significant RNA marker of some cluster, AND does its ATAC gene-activity "
            "independently confirm elevated accessibility in that cluster's real cross-modal "
            "partner (matched by cell overlap)? Use this to check whether a specific marker "
            "you're relying on (from your own analysis or from literature) is corroborated by "
            "BOTH modalities independently, not just RNA."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"gene": {"type": "string", "description": "Gene symbol"}},
            "required": ["gene"],
        },
    },
    {
        "name": "check_for_condition_groups",
        "description": (
            "Check whether this dataset has a drug/treatment/condition axis (e.g. drug vs. "
            "vehicle control) beyond cell type/cell line -- by column name, not by guessing. "
            "Call this once, early, before assuming every dataset is control-only. If it finds "
            "one, follow up with `condition_group_qc` before analyzing conditions separately."
        ),
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "condition_group_qc",
        "description": (
            "Per-arm cell counts for one condition column found by check_for_condition_groups "
            "-- flags any arm with too few cells to support a condition-level claim."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "column_key": {"type": "string", "description": "e.g. 'rna.drug_treatment', from check_for_condition_groups' output"},
            },
            "required": ["column_key"],
        },
    },
    {
        "name": "get_qc_summary",
        "description": (
            "Get a summary of this dataset's fixed-core QC and clustering numbers "
            "(genes/UMIs/mito/doublets per cell, ATAC fragments/FRiP/TSS enrichment/"
            "nucleosome signal, cluster counts, RNA-ATAC cluster agreement). Call this "
            "first when asked to assess data quality -- don't guess at QC numbers, look "
            "them up."
        ),
        "input_schema": {"type": "object", "properties": {}},
    },
]

# A deliberate, bounded exception to the MVP's "agent does not write new
# analysis code" rule (CLAUDE.md) -- kept OUT of `TOOLS` above so it's opt-in
# per investigation (`tools=TOOLS_WITH_CODE_EXECUTION`), not silently
# available to every run. Only ever backed by a small, pre-aggregated
# `code_namespace` dict the caller builds ahead of time -- never the raw
# per-cell mdata -- see `tools/code_execution.py` for the sandboxing.
CODE_EXECUTION_TOOL = {
    "name": "run_analysis_code",
    "description": (
        "Run a short Python snippet (pandas as `pd`, numpy as `np`) over a small "
        "pre-aggregated summary dict available as `data` -- NOT the raw per-cell "
        "dataset, only whatever aggregate statistics were provided for this "
        "investigation. Assign your answer to a variable named `result`, or use "
        "print(). No imports, no file/network/system access -- these are rejected."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "code": {"type": "string", "description": "Python snippet using `pd`, `np`, `data`"},
        },
        "required": ["code"],
    },
}

TOOLS_WITH_CODE_EXECUTION = TOOLS + [CODE_EXECUTION_TOOL]


@dataclass
class AgentRunResult:
    question: str
    answer: str
    turn_count: int
    hit_max_turns: bool
    input_tokens: int
    output_tokens: int
    cache_read_tokens: int
    cache_creation_tokens: int
    estimated_cost_usd: float
    reasoning_trail_path: str
    tool_calls: list = field(default_factory=list)
    refused: bool = False
    refusal_detail: str | None = None


class ReasoningTrail:
    """Writes a human-readable markdown log of the run as it happens."""

    def __init__(self, question: str, model: str):
        RUNS_DIR.mkdir(parents=True, exist_ok=True)
        run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        self.path = RUNS_DIR / f"{run_id}.md"
        self.path.write_text(f"# Agent run {run_id}\n\n**Model:** {model}\n\n**Question:** {question}\n\n---\n\n")

    def log_turn(self, turn: int, text_blocks: list[str], tool_calls: list[dict]):
        with self.path.open("a") as f:
            f.write(f"## Turn {turn}\n\n")
            for text in text_blocks:
                f.write(f"{text}\n\n")
            for call in tool_calls:
                f.write(f"**Tool call:** `{call['name']}({json.dumps(call['input'])})`\n\n")
                f.write(f"**Result:** {call['result_summary']}\n\n")

    def log_final(self, answer: str, hit_max_turns: bool, cost: float, turns: int):
        with self.path.open("a") as f:
            f.write("---\n\n## Final answer\n\n")
            f.write(answer + "\n\n")
            if hit_max_turns:
                f.write(f"*(Stopped: hit max-turns limit at turn {turns}.)*\n\n")
            f.write(f"**Turns:** {turns} | **Estimated cost:** ${cost:.5f}\n")


_RESULT_SUMMARY_MAX_CHARS = 1000  # the QC summary tool's real output alone
# already runs ~340 chars (500 cells, 15695 genes, ... ARI: 0.467.) even
# before an eval scenario's fault-specific note is appended -- 300 truncated
# it mid-sentence, which would've silently dropped real numbers (like the
# ARI figure) from what the grounding checker can see, causing spurious
# "ungrounded" flags on citations that actually did come from a real tool
# result. Caught by measuring the real string length before running the
# eval harness for real, not by assumption.


def _summarize_result(result) -> str:
    text = json.dumps(result) if not isinstance(result, str) else result
    if len(text) <= _RESULT_SUMMARY_MAX_CHARS:
        return text
    return text[:_RESULT_SUMMARY_MAX_CHARS] + "... (truncated for log)"


def _execute_tool(name: str, tool_input: dict, mdata, qc_summary: str | None, code_namespace: dict | None = None) -> tuple[str, bool]:
    """Returns (result_content_str, is_error). `qc_summary` is computed
    lazily here (only when the tool is actually invoked), not eagerly in
    `run_agent` -- eagerly calling `fixed_core_summary(mdata)` would crash
    whenever `mdata` is None (e.g. the max-turns-guard test, which never
    exercises real tools and monkeypatches `get_agent_fixed_core` to return
    None), even on runs that never call this tool at all."""
    try:
        if name == "tf_motif_correlation":
            result = tf_motif_correlation(mdata, tool_input["gene"])
        elif name == "peak_to_gene_links":
            result = peak_to_gene_links(mdata, tool_input["gene"])
        elif name == "regulon_inference":
            result = regulon_inference(mdata, tool_input["tf_gene"], target_gene=tool_input.get("target_gene"))
        elif name == "search_pubmed":
            result = search_pubmed(tool_input["query"])
        elif name == "fetch_pubmed_abstracts":
            result = fetch_pubmed_abstracts(tool_input["pmids"])
        elif name == "enrich_gene_set":
            result = enrich_gene_set(tool_input["genes"])
        elif name == "resolve_depmap_id":
            result = resolve_depmap_id(tool_input["depmap_id"])
        elif name == "check_for_identity_columns":
            result = detect_identity_columns(mdata)
        elif name == "list_clusters":
            result = {
                "rna_clusters": sorted(mdata.mod["rna"].obs[RNA_CLUSTER_KEY].unique().tolist(), key=str),
                "atac_clusters": sorted(mdata.mod["atac"].obs[ATAC_CLUSTER_KEY].unique().tolist(), key=str),
            }
        elif name == "top_cluster_markers":
            result = top_rna_cluster_markers(mdata.mod["rna"], tool_input["cluster"], n=tool_input.get("n", 15))
        elif name == "top_gene_activity_markers":
            result = top_gene_activity_cluster_markers(mdata.mod["atac"], tool_input["cluster"], n=tool_input.get("n", 15))
        elif name == "cross_modal_marker_check":
            result = cross_modal_marker_validation(mdata, tool_input["gene"])
        elif name == "record_checklist_item":
            return json.dumps(tool_input), False  # echoed back; the caller extracts structured items from tool_calls
        elif name == "record_novel_finding":
            return json.dumps(tool_input), False
        elif name == "record_judger_verdict":
            return json.dumps(tool_input), False
        elif name == "check_for_condition_groups":
            result = detect_condition_groups(mdata)
        elif name == "condition_group_qc":
            result = condition_group_qc(mdata, tool_input["column_key"])
        elif name == "get_qc_summary":
            return qc_summary if qc_summary is not None else format_qc_summary(fixed_core_summary(mdata)), False
        elif name == "run_analysis_code":
            return run_analysis_code(tool_input["code"], code_namespace or {}), False
        else:
            return f"Unknown tool: {name}", True
        return json.dumps(result), False
    except Exception as e:
        logger.exception("Tool %s failed", name)
        return f"Tool error: {e}", True


def run_agent(
    question: str,
    model: str | None = None,
    max_turns: int | None = None,
    mdata=None,
    qc_summary: str | None = None,
    tools: list | None = None,
    code_namespace: dict | None = None,
    dataset_context: str | None = None,
    require_tool_call: str | None = None,
) -> AgentRunResult:
    """`mdata` and `qc_summary` default to the cached clean fixed-core result
    and its own real summary (unchanged behavior for existing callers). The
    eval harness passes both explicitly per fault scenario, so `get_qc_summary`
    reflects that scenario's real (possibly faulted) numbers rather than
    always the clean baseline -- letting the agent's claims stay traceable to
    a real tool result (and thus grounding-checkable) instead of numbers
    baked into the prompt where they'd bypass the grounding check entirely.

    `tools` defaults to the standard 4-tool `TOOLS` (unchanged behavior for
    existing callers) -- pass `TOOLS_WITH_CODE_EXECUTION` explicitly to opt a
    specific investigation into the sandboxed code-execution tool, paired
    with a `code_namespace` (pre-aggregated statistics only, never raw
    per-cell data).

    `dataset_context` defaults to `PBMC_DATASET_CONTEXT` (unchanged behavior
    for every pre-existing caller) -- pass `OWN_DATA_CONTEXT` (from
    `agent.prompts`) for shareseq-multi-cell-lines runs. Fixes a real bug: before this
    parameter existed, every run used the PBMC-specific system-prompt text
    unconditionally, and a real eval run against the shareseq-multi-cell-lines data caught a
    model (Opus) explicitly flagging the resulting "premise mismatch" (told
    it's PBMC, tool results say otherwise) as a data-quality problem in its
    own right -- not a hypothetical concern, an actual observed failure.

    `require_tool_call`, if set, names a tool that MUST be called before the
    run is allowed to end on `end_turn`. Real bug this fixes: the Judger's
    TCF7L2-motif-degeneracy verdict (a genuinely hard case) spiralled through
    ~15 `search_pubmed` queries, ran out of promising leads, and ended with a
    long prose conclusion instead of calling `record_judger_verdict` -- the
    whole finding silently rendered as "NO VERDICT RECORDED" in the report.
    One forced-`tool_choice` nudge turn recovers it instead of losing the
    verdict outright.
    """
    model = model or AGENT_MODEL
    max_turns = max_turns or MAX_TURNS
    tools = tools if tools is not None else TOOLS
    dataset_context = dataset_context if dataset_context is not None else PBMC_DATASET_CONTEXT

    client = anthropic.Anthropic(api_key=get_anthropic_api_key())
    mdata = mdata if mdata is not None else get_agent_fixed_core()
    trail = ReasoningTrail(question, model)

    messages = [{"role": "user", "content": question}]
    total_in = total_out = total_cache_read = total_cache_creation = 0
    all_tool_calls = []
    turn = 0
    hit_max_turns = False
    final_text = ""
    refused = False
    refusal_detail = None
    force_tool = None  # set to require_tool_call for exactly one nudge turn
    nudged = False

    logger.info("Starting agent run: model=%s max_turns=%d question=%r", model, max_turns, question)

    while True:
        if turn >= max_turns:
            hit_max_turns = True
            logger.warning("Hit max_turns=%d, stopping loop", max_turns)
            final_text = final_text or "(Stopped: reached max-turns limit before producing a final answer.)"
            break
        turn += 1

        try:
            response = client.messages.create(
                model=model,
                # 1024 was fine for Haiku's terse style but silently truncated
                # Sonnet/Opus/Fable mid-answer on 4/12 real eval runs (verified:
                # empty answer text, stop_reason="max_tokens") -- these models
                # write more elaborate multi-paragraph answers by default even
                # for the same question. 4096 gives real headroom without
                # inflating cost much (this loop's answers run a few hundred
                # tokens even when NOT truncated).
                max_tokens=4096,
                system=[
                    {"type": "text", "text": CORE_SYSTEM_PROMPT, "cache_control": {"type": "ephemeral"}},
                    {"type": "text", "text": dataset_context},
                ],
                tools=tools,
                messages=messages,
                **({"tool_choice": {"type": "tool", "name": force_tool}} if force_tool else {}),
            )
            force_tool = None
        except anthropic.RateLimitError as e:
            logger.error("Rate limited: %s", e)
            raise
        except anthropic.APIStatusError as e:
            logger.error("API status error: %s", e)
            raise
        except anthropic.APIConnectionError as e:
            logger.error("API connection error: %s", e)
            raise

        usage = response.usage
        total_in += usage.input_tokens
        total_out += usage.output_tokens
        total_cache_read += getattr(usage, "cache_read_input_tokens", 0) or 0
        total_cache_creation += getattr(usage, "cache_creation_input_tokens", 0) or 0

        text_blocks = [b.text for b in response.content if b.type == "text"]
        tool_use_blocks = [b for b in response.content if b.type == "tool_use"]

        if response.stop_reason == "end_turn":
            final_text = "\n".join(text_blocks) or final_text
            trail.log_turn(turn, text_blocks, [])
            already_called = any(
                tc["name"] == require_tool_call and not tc.get("is_error", False) for tc in all_tool_calls
            )
            if require_tool_call and not already_called and not nudged:
                nudged = True
                force_tool = require_tool_call
                messages.append({"role": "assistant", "content": response.content})
                messages.append({"role": "user", "content": (
                    f"You stopped without calling `{require_tool_call}`. Based on everything "
                    f"you've already found, call `{require_tool_call}` now with your final answer."
                )})
                logger.warning("Agent ended without required tool %s at turn %d, nudging once", require_tool_call, turn)
                continue
            logger.info("Agent finished at turn %d", turn)
            break

        if response.stop_reason == "refusal":
            # Track honestly rather than retry-until-compliant or silently
            # drop -- CLAUDE.md wants refusals/hedging visible as their own
            # signal in the model-comparison table, especially since Fable 5
            # applies extra safety classifiers to biology content that other
            # models in this comparison don't. `stop_details` is populated
            # only on a refusal (null otherwise) per the API's own contract.
            details = getattr(response, "stop_details", None)
            category = getattr(details, "category", None) if details else None
            explanation = getattr(details, "explanation", None) if details else None
            refused = True
            refusal_detail = f"category={category!r} explanation={explanation!r}"
            final_text = "\n".join(text_blocks) or f"(Model declined to answer: {refusal_detail})"
            trail.log_turn(turn, text_blocks or [f"[REFUSAL: {refusal_detail}]"], [])
            logger.warning("Agent run refused at turn %d: %s", turn, refusal_detail)
            break

        messages.append({"role": "assistant", "content": response.content})

        turn_calls = []
        tool_results = []
        for block in tool_use_blocks:
            result_str, is_error = _execute_tool(block.name, block.input, mdata, qc_summary, code_namespace)
            logger.info("Tool call: %s(%s) -> error=%s", block.name, block.input, is_error)
            tool_results.append({
                "type": "tool_result",
                "tool_use_id": block.id,
                "content": result_str,
                "is_error": is_error,
            })
            turn_calls.append({
                "name": block.name, "input": block.input,
                # "result" is the FULL, untruncated tool output -- grounding
                # checks need this (see eval/grounding.py); a fetch_pubmed_abstracts
                # call alone can return several KB across multiple abstracts,
                # well past _RESULT_SUMMARY_MAX_CHARS, so relying on the
                # truncated summary here would silently hide later PMIDs/
                # numbers from the grounding checker (the exact failure mode
                # that made _RESULT_SUMMARY_MAX_CHARS get raised once before
                # -- this time the fix is to stop truncating the field
                # grounding reads at all, not raise the cap again).
                # "result_summary" stays truncated -- only the human-readable
                # reasoning-trail log needs a short display string.
                "result": result_str, "result_summary": _summarize_result(result_str), "is_error": is_error,
            })
        all_tool_calls.extend(turn_calls)
        trail.log_turn(turn, text_blocks, turn_calls)

        if tool_results:
            messages.append({"role": "user", "content": tool_results})
        elif response.stop_reason != "tool_use":
            # No tool calls and not end_turn (e.g. max_tokens) -- stop rather than loop forever.
            final_text = "\n".join(text_blocks) or f"(Stopped: unexpected stop_reason={response.stop_reason!r}.)"
            break

    rate = MODEL_COSTS.get(model, {"input": 0.0, "output": 0.0})
    # `usage.input_tokens` is already exclusive of cached tokens (Anthropic
    # reports cache_creation/cache_read as separate, non-overlapping counts),
    # so it's charged at the full input rate with no subtraction here.
    cost = total_in / 1e6 * rate["input"] + total_out / 1e6 * rate["output"]
    cost += total_cache_creation / 1e6 * rate["input"] * 1.25
    cost += total_cache_read / 1e6 * rate["input"] * 0.1

    trail.log_final(final_text, hit_max_turns, cost, turn)
    logger.info(
        "Agent run complete: turns=%d cost=$%.5f in=%d out=%d cache_read=%d",
        turn, cost, total_in, total_out, total_cache_read,
    )

    return AgentRunResult(
        question=question, answer=final_text, turn_count=turn, hit_max_turns=hit_max_turns,
        input_tokens=total_in, output_tokens=total_out,
        cache_read_tokens=total_cache_read, cache_creation_tokens=total_cache_creation,
        estimated_cost_usd=cost, reasoning_trail_path=str(trail.path), tool_calls=all_tool_calls,
        refused=refused, refusal_detail=refusal_detail,
    )


if __name__ == "__main__":
    q = " ".join(sys.argv[1:]) or "What can you tell me about SPI1 in this dataset?"
    result = run_agent(q)
    print(result.answer)
    print(f"\n[turns={result.turn_count} cost=${result.estimated_cost_usd:.5f} trail={result.reasoning_trail_path}]")
