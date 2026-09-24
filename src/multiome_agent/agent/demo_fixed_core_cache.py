"""A small, deliberately-trimmed fixed-core object committed directly to
the repo (`demo_cache/agent_demo_core.h5mu`, ~46MB) so `make demo` is
genuinely fast on a fresh clone -- unlike the full report-generation
pipeline (`agent/fixed_core_cache.py`'s `get_agent_fixed_core`: 11,909
cells, ~2.6GB, never committed, since GitHub hard-rejects any git-tracked
file over 100MB and this project's own real full-scale rebuild takes
60-90+ minutes).

Committable at all only because it's trimmed to exactly what the 11 agent
tools (`agent/loop.py`'s `TOOLS`) read at query time -- confirmed by
grepping every tool's implementation, not guessed: RNA log-normalized
expression + QC/cluster `.obs` columns + `rank_genes_groups`; ATAC
QC/cluster `.obs` columns + gene-activity + chromVAR-deviations matrices.
Dropped: PCA/LSI loadings and descriptive peak/gene annotation columns
(clustering scratch needed only to BUILD the cluster labels/DE tables that
are kept, never read again afterward) and the raw ATAC peak matrix (also
unused by any tool once gene activity/chromVAR are computed from it). This
shrunk a naive n=150 rebuild from 96.8MB to 46.25MB. Also strips
`atac.uns["files"]`, which otherwise bakes this machine's absolute local
fragments-file path into the committed binary -- a real leak caught before
committing, not a hypothetical one.

Illustrative only, at a small subsample deliberately not meant to support
rigorous claims -- see `reports/examples/` for the real, full-scale
analyses this project's flagship reports are built from.
"""

from __future__ import annotations

import scipy.sparse as sp
from mudata import MuData, read_h5mu, write_h5mu

from multiome_agent.config import REPO_ROOT
from multiome_agent.core.pipeline import run_fixed_core
from multiome_agent.data.loader import RAW_DIR, load_pbmc_multiome
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

DEMO_N_CELLS = 150
DEMO_SEED = 0
FRAGMENTS_FILE = RAW_DIR / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
DEMO_CACHE_PATH = REPO_ROOT / "demo_cache" / "agent_demo_core.h5mu"


def _trim_to_tool_readable_fields(mdata: MuData) -> MuData:
    """Keeps only what `agent/loop.py`'s `TOOLS` actually read at query
    time (see module docstring) -- everything else is clustering/embedding
    scratch or purely descriptive annotation that's dead weight once
    cluster labels and DE/gene-activity/chromVAR results are already baked
    in, plus one real path leak (`atac.uns["files"]`)."""
    rna, atac = mdata.mod["rna"], mdata.mod["atac"]

    rna.var = rna.var[[]]
    rna.varm.clear()
    rna.obsm.clear()
    for key in list(rna.uns.keys()):
        if key != "rank_genes_groups":
            del rna.uns[key]
    rna.layers.pop("counts", None)

    atac.var = atac.var[[]]
    atac.varm.clear()
    for key in ("X_lsi", "X_lsi_no_depth"):
        atac.obsm.pop(key, None)
    for key in list(atac.uns.keys()):
        if key not in ("gene_activity_genes", "chromvar_motif_names"):
            del atac.uns[key]  # drops uns["files"]'s local fragments path too
    atac.layers.pop("counts", None)
    atac.X = sp.csr_matrix((atac.n_obs, atac.n_vars))  # unused by any tool; kept shape-valid, zero storage cost

    mdata.update()
    return mdata


def build_demo_cache(n_cells: int = DEMO_N_CELLS, seed: int = DEMO_SEED) -> MuData:
    """Builds the small, trimmed, committable demo cache from scratch and
    writes it to `DEMO_CACHE_PATH`. Real work -- still pays the ~17-20 min
    fixed fragments-import cost (see `core/pipeline.py`) -- run once by a
    maintainer and committed, never run by end users. See
    `scripts/build_demo_cache.py` for the CLI entry point.
    """
    mdata = load_pbmc_multiome(n_cells=n_cells, seed=seed, force_reload=False)
    mdata = run_fixed_core(mdata, FRAGMENTS_FILE)
    mdata = _trim_to_tool_readable_fields(mdata)
    DEMO_CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    write_h5mu(DEMO_CACHE_PATH, mdata)
    logger.info("Wrote demo cache to %s", DEMO_CACHE_PATH)
    return mdata


def get_demo_fixed_core() -> MuData:
    """Loads the committed demo cache. Errors clearly if missing rather
    than silently rebuilding -- rebuilding would defeat the entire point
    (`make demo` must never trigger the ~17-20 min pipeline on an end
    user's machine)."""
    if not DEMO_CACHE_PATH.exists():
        raise FileNotFoundError(
            f"{DEMO_CACHE_PATH} not found. This is meant to be committed to the repo -- "
            "if you're intentionally regenerating it, run scripts/build_demo_cache.py."
        )
    return read_h5mu(DEMO_CACHE_PATH)
