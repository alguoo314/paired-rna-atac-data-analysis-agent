"""Disk cache for the shareseq-multi-cell-lines fixed-core result (now
including chromVAR motif deviations, added after the initial Stretch-phase
build wrongly deferred them -- see PROGRESS.md's correction entry).

The full run (5,814 cells, ~134K filtered peaks x 879 JASPAR motifs) takes
~39 minutes, dominated by chromVAR's background-deviation permutation step --
mirrors the tenx-cell-ranger pipeline's `agent/fixed_core_cache.py` pattern
(build once, reuse via `.h5mu` on disk) for exactly the same reason: too slow
to recompute per tool call or test run.
"""

from __future__ import annotations

from pathlib import Path

from mudata import MuData, read_h5mu, write_h5mu

from multiome_agent.config import REPO_ROOT
from multiome_agent.core.shareseq_pipeline import run_shareseq_fixed_core
from multiome_agent.data.shareseq_loader import load_shareseq_multiome
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

CACHE_PATH = REPO_ROOT / "data" / "processed" / "shareseq_fixed_core.h5mu"


def get_shareseq_fixed_core(force_reload: bool = False) -> MuData:
    """Load-or-build the shareseq-multi-cell-lines data's fixed-core result
    (QC, clustering, cross-modal agreement, cell-line recovery, motif
    deviations)."""
    if CACHE_PATH.exists() and not force_reload:
        logger.info("Shareseq fixed-core cache hit: %s", CACHE_PATH)
        return read_h5mu(CACHE_PATH)

    logger.info("Shareseq fixed-core cache miss: building (~39 min, dominated by chromVAR)")
    mdata = load_shareseq_multiome()
    mdata = run_shareseq_fixed_core(mdata)

    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    write_h5mu(CACHE_PATH, mdata)
    logger.info("Wrote shareseq fixed-core cache to %s", CACHE_PATH)
    return mdata
