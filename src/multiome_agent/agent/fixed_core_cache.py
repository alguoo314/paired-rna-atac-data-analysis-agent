"""Fixed-core result the agent's tools operate on, and what the
comprehensive report is built from.

Previously subsampled to n_cells=500 for iteration speed -- but that made
the committed report's own fault-injection eval unreliable: models
correctly noticed 500 cells split across 8-11 Leiden clusters is too few
per cluster and repeatedly flagged it, even on the clean control, as a real
(if second-order) limitation rather than the injected fault. Fixed by using
the full dataset (11,909 cells, the real count from the raw Cell Ranger ARC
output -- see `data/loader.py`'s inspection notes) instead of a subsample.
This does NOT multiply the dominant cost: snapatac2's fragments-file
sort/import is a ~17-minute cost fixed by the fragments file's size, not by
how many cells are kept (see PROGRESS.md) -- gene-activity/chromVAR scale
with cell count, but from a substantially smaller base than the fixed cost.
Still slow enough to disk-cache once and reuse rather than recompute per
tool call or test run.
"""

from __future__ import annotations

from pathlib import Path

from mudata import MuData, read_h5mu, write_h5mu

from multiome_agent.config import REPO_ROOT
from multiome_agent.core.pipeline import run_fixed_core
from multiome_agent.data.loader import RAW_DIR, load_pbmc_multiome
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

AGENT_N_CELLS = 11909  # the real, full cell count -- not a subsample
AGENT_SEED = 0
FRAGMENTS_FILE = RAW_DIR / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
CACHE_PATH = REPO_ROOT / "data" / "processed" / f"agent_fixed_core_n{AGENT_N_CELLS}_seed{AGENT_SEED}.h5mu"


def get_agent_fixed_core(force_reload: bool = False) -> MuData:
    """Load-or-build the small fixed-core result the agent's tools use."""
    if CACHE_PATH.exists() and not force_reload:
        logger.info("Agent fixed-core cache hit: %s", CACHE_PATH)
        return read_h5mu(CACHE_PATH)

    logger.info("Agent fixed-core cache miss: building (n_cells=%d)", AGENT_N_CELLS)
    mdata = load_pbmc_multiome(n_cells=AGENT_N_CELLS, seed=AGENT_SEED, force_reload=False)
    mdata = run_fixed_core(mdata, FRAGMENTS_FILE)

    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    write_h5mu(CACHE_PATH, mdata)
    logger.info("Wrote agent fixed-core cache to %s", CACHE_PATH)
    return mdata
