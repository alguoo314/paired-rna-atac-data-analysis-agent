"""Shared fixtures for fixed-core pipeline tests.

The full fixed-core run (RNA/ATAC QC, clustering, gene activity, motif
deviations) is expensive -- dominated by a ~3.5 min snapatac2 fragments-file
import -- so it runs at most once per test session via this session-scoped
fixture, shared across test_core_analyses.py and
test_gene_activity_and_motifs.py.
"""

import pytest

from multiome_agent.core.pipeline import run_fixed_core
from multiome_agent.data.loader import (
    MATRIX_H5,
    PEAK_ANNOTATION_TSV,
    PER_BARCODE_METRICS_CSV,
    RAW_DIR,
    load_pbmc_multiome,
)

FRAGMENTS_FILE = RAW_DIR / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
GENOME_2BIT = RAW_DIR.parent / "reference" / "hg38.2bit"

RAW_DATA_PRESENT = (
    MATRIX_H5.exists()
    and PER_BARCODE_METRICS_CSV.exists()
    and PEAK_ANNOTATION_TSV.exists()
    and FRAGMENTS_FILE.exists()
)
FULL_FIXED_CORE_DATA_PRESENT = RAW_DATA_PRESENT and GENOME_2BIT.exists()

N_CELLS = 3000


@pytest.fixture(scope="session")
def fixed_core_mdata():
    if not FULL_FIXED_CORE_DATA_PRESENT:
        pytest.skip(
            "raw PBMC multiome files (incl. ATAC fragments + hg38.2bit) "
            "not present under data/raw/"
        )
    mdata = load_pbmc_multiome(n_cells=N_CELLS, seed=0, force_reload=False)
    return run_fixed_core(mdata, FRAGMENTS_FILE)


@pytest.fixture(scope="session")
def agent_fixed_core_mdata():
    """The agent's smaller, faster fixed-core cache (see agent/fixed_core_cache.py) --
    separate from `fixed_core_mdata` above (n=3000, eval-quality) since agent tool
    tests don't need eval-scale data, just a real, already-annotated MuData."""
    if not FULL_FIXED_CORE_DATA_PRESENT:
        pytest.skip(
            "raw PBMC multiome files (incl. ATAC fragments + hg38.2bit) "
            "not present under data/raw/"
        )
    from multiome_agent.agent.fixed_core_cache import get_agent_fixed_core

    return get_agent_fixed_core()
