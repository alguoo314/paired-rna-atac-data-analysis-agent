"""Small, fast fixed-core result the agent's tools operate on.

The eval-quality fixed-core run uses n_cells=3000 and takes ~17 minutes
(dominated by snapatac2's fragments-file import + motif scanning) -- too
slow to recompute per tool call or test run. This module builds a much
smaller subsample once (n_cells=500, still large enough for clusters/
markers/motifs to be meaningful) and caches the full fixed-core result
(QC, clustering, gene activity, motif deviations) as .h5mu.
"""

from __future__ import annotations

from pathlib import Path

from mudata import MuData, read_h5mu, write_h5mu

from multiome_agent.config import REPO_ROOT
from multiome_agent.core.pipeline import run_fixed_core
from multiome_agent.data.loader import RAW_DIR, load_pbmc_multiome
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

AGENT_N_CELLS = 500
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
