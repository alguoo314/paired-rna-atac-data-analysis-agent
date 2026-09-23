"""Maps an agent-decided loading strategy (see `agent.loader_selector`) to
the actual deterministic loader function. The agent only ever picks between
these two named strategies via a schema-constrained tool call -- it never
writes or chooses parsing code itself.
"""

from __future__ import annotations

from mudata import MuData

from multiome_agent.agent.loader_selector import LoaderDecision, decide_loading_strategy
from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD
from multiome_agent.data.loader import MATRIX_H5, load_pbmc_multiome
from multiome_agent.data.shareseq_loader import load_shareseq_multiome
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def candidate_paths_for(source: str) -> list[str]:
    """`source` is an operational choice by whoever is running a report
    (which data location to point at this time), NOT a loading-strategy
    choice -- that part is still left to the agent. Raw file paths only, so
    the agent has real structure to inspect rather than a name to trust.
    """
    if source == "tenx-cell-ranger":
        return [str(MATRIX_H5)]
    if source == "shareseq-multi-cell-lines":
        return [SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD, SHARESEQ_ATAC_H5AD]
    raise ValueError(f"Unknown source {source!r}, expected 'tenx-cell-ranger' or 'shareseq-multi-cell-lines'")


def load_dataset_via_agent_decision(source: str, model: str | None = None) -> tuple[MuData, LoaderDecision]:
    """Full agent-mediated load: inspect -> decide -> deterministically load.
    `source` only selects WHICH candidate paths to hand the agent -- the
    agent still independently determines the loading strategy from their
    real structure (see `agent.loader_selector.decide_loading_strategy`).
    """
    paths = candidate_paths_for(source)
    decision = decide_loading_strategy(paths, model=model)

    if decision.strategy == "combined_single_file":
        mdata = load_pbmc_multiome()
    elif decision.strategy == "separate_per_modality_files":
        mdata = load_shareseq_multiome()
    else:
        raise ValueError(f"Unknown loading strategy decided: {decision.strategy!r}")

    logger.info(
        "Loaded dataset for source=%r via agent-decided strategy=%s (%d cells)",
        source, decision.strategy, mdata.n_obs,
    )
    return mdata, decision


def load_fixed_core_via_agent_decision(source: str, model: str | None = None) -> tuple[MuData, LoaderDecision]:
    """Same agent-mediated loading-strategy decision as
    `load_dataset_via_agent_decision`, but hands off to the CACHED,
    already fixed-core-processed getters (gene activity + chromVAR already
    computed, on-disk cache) rather than the raw loaders -- what the report
    generator actually needs, without repaying the expensive fixed-core
    computation (chromVAR alone: ~39 min for the shareseq-multi-cell-lines data) that's
    already been paid for and cached on disk.
    """
    paths = candidate_paths_for(source)
    decision = decide_loading_strategy(paths, model=model)

    if decision.strategy == "combined_single_file":
        from multiome_agent.agent.fixed_core_cache import get_agent_fixed_core
        mdata = get_agent_fixed_core()
    elif decision.strategy == "separate_per_modality_files":
        from multiome_agent.agent.shareseq_fixed_core_cache import get_shareseq_fixed_core
        mdata = get_shareseq_fixed_core()
    else:
        raise ValueError(f"Unknown loading strategy decided: {decision.strategy!r}")

    logger.info(
        "Loaded CACHED fixed-core result for source=%r via agent-decided strategy=%s (%d cells)",
        source, decision.strategy, mdata.n_obs,
    )
    return mdata, decision
