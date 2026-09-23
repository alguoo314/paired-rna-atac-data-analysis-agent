"""Agent-mediated loader selection: given a list of raw data file paths, a
small dedicated agent run inspects their real structure and decides whether
this is the "combined single file" format (tenx-cell-ranger-style, one file with
both RNA and ATAC features) or the "separate per-modality files" format
(shareseq-style, one file per modality). The decision is a required,
schema-constrained tool call, not free text -- so parsing it back out is
exact, not a regex over prose.

Deliberately a separate, minimal loop rather than reusing `agent.loop.run_agent`:
that loop's `CORE_SYSTEM_PROMPT` narrates PBMC-flavored analysis tools
(`tf_motif_correlation`, `search_pubmed`, ...) that have nothing to do with
this bootstrap step, and threading a `require_fixed_core`-style escape hatch
through its "PBMC fixed-core cache" default would add more coupling than
this two-tool, few-turn task is worth. Once the format is decided, the
ACTUAL parsing is fully deterministic Python (`load_pbmc_multiome` /
`load_shareseq_multiome`) -- only the choice is agent-mediated, per
CLAUDE.md's "LLM never touches raw data" principle.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import anthropic

from multiome_agent.config import AGENT_MODEL, get_anthropic_api_key
from multiome_agent.logging_utils import get_logger
from multiome_agent.tools.data_inspection import inspect_data_paths

logger = get_logger(__name__)

# Mirrors agent/loop.py's MODEL_COSTS -- duplicated rather than imported to
# keep this module fully independent of agent.loop (see module docstring);
# both must be updated together if pricing changes.
MODEL_COSTS = {
    "claude-haiku-4-5": {"input": 1.00, "output": 5.00},
    "claude-sonnet-5": {"input": 3.00, "output": 15.00},
    "claude-opus-5": {"input": 5.00, "output": 25.00},
    "claude-fable-5": {"input": 10.00, "output": 50.00},
}

SYSTEM_PROMPT = """\
You determine how a set of raw single-cell multiome data files should be \
loaded, based only on their real on-disk structure -- never assume the \
format from a file name or path.

Two loading strategies exist:
- "combined_single_file": ONE file contains both RNA and ATAC features \
together (its feature-type list includes both "Gene Expression" and "Peaks").
- "separate_per_modality_files": RNA and ATAC live in DIFFERENT files, each \
holding only one modality.

Call `inspect_data_paths` with the given file paths to see their real \
structure (format, shape, feature types, column names). Then call \
`decide_loading_strategy` exactly once with your decision and the specific \
structural evidence that led to it.
"""

INSPECT_TOOL = {
    "name": "inspect_data_paths",
    "description": (
        "Peek at the structure of one or more raw data files (10x Cell Ranger "
        ".h5 or AnnData .h5ad) without loading their full contents. Returns, "
        "per path: file format, number of observations/variables, column "
        "names, and (for combined-format files) which feature types are "
        "present."
    ),
    "input_schema": {
        "type": "object",
        "properties": {"paths": {"type": "array", "items": {"type": "string"}}},
        "required": ["paths"],
    },
}

DECIDE_TOOL = {
    "name": "decide_loading_strategy",
    "description": "Report your final decision on how to load this data source.",
    "input_schema": {
        "type": "object",
        "properties": {
            "strategy": {
                "type": "string",
                "enum": ["combined_single_file", "separate_per_modality_files"],
            },
            "reasoning": {
                "type": "string",
                "description": "One or two sentences citing the specific structural evidence (e.g. feature types seen, number of files, shapes) that led to this decision.",
            },
        },
        "required": ["strategy", "reasoning"],
    },
}


@dataclass
class LoaderDecision:
    strategy: str
    reasoning: str
    cost_usd: float
    turns: int


def _execute(name: str, tool_input: dict) -> tuple[str, bool]:
    try:
        if name == "inspect_data_paths":
            return json.dumps(inspect_data_paths(tool_input["paths"])), False
        if name == "decide_loading_strategy":
            return json.dumps(tool_input), False
        return f"Unknown tool: {name}", True
    except Exception as e:
        logger.exception("Tool %s failed", name)
        return f"Tool error: {e}", True


def decide_loading_strategy(paths: list[str], model: str | None = None, max_turns: int = 4) -> LoaderDecision:
    """Runs the bootstrap decision loop and returns the parsed decision --
    does NOT load the data itself (see `load_dataset_via_agent_decision`,
    which does, using this function's result).
    """
    model = model or AGENT_MODEL
    client = anthropic.Anthropic(api_key=get_anthropic_api_key())
    messages = [{"role": "user", "content": f"File paths for this run: {json.dumps(paths)}"}]
    total_in = total_out = 0
    decision = None

    for turn in range(1, max_turns + 1):
        response = client.messages.create(
            model=model,
            max_tokens=1024,
            system=SYSTEM_PROMPT,
            tools=[INSPECT_TOOL, DECIDE_TOOL],
            messages=messages,
        )
        total_in += response.usage.input_tokens
        total_out += response.usage.output_tokens

        tool_use_blocks = [b for b in response.content if b.type == "tool_use"]
        if response.stop_reason != "tool_use" or not tool_use_blocks:
            break

        messages.append({"role": "assistant", "content": response.content})
        tool_results = []
        for block in tool_use_blocks:
            result_str, is_error = _execute(block.name, block.input)
            tool_results.append({
                "type": "tool_result", "tool_use_id": block.id,
                "content": result_str, "is_error": is_error,
            })
            if block.name == "decide_loading_strategy" and not is_error:
                decision = block.input
        messages.append({"role": "user", "content": tool_results})

        if decision is not None:
            break

    if decision is None:
        raise RuntimeError(
            f"Loader-selection agent never called decide_loading_strategy within {max_turns} turns "
            f"for paths={paths!r} -- can't determine loading strategy."
        )

    rate = MODEL_COSTS.get(model, {"input": 0.0, "output": 0.0})
    cost = total_in / 1e6 * rate["input"] + total_out / 1e6 * rate["output"]
    logger.info(
        "Loader decision for %r: strategy=%s cost=$%.5f (%s)",
        paths, decision["strategy"], cost, decision["reasoning"],
    )
    return LoaderDecision(strategy=decision["strategy"], reasoning=decision["reasoning"], cost_usd=cost, turns=turn)
