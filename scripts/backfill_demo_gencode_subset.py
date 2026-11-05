"""One-time repair for the already-committed demo cache
(`demo_cache/agent_demo_core.h5mu`), backfilling
`mdata.uns["gencode_gene_coords_subset"]` -- added after the cache was
first built (see `agent/demo_fixed_core_cache.py`'s `_embed_gencode_gene_
coords_subset`, which a fresh `build_demo_cache()` run would now do on its
own). Deliberately NOT a rebuild from raw files (that pays the ~17-20 min
fragments-reimport cost again for data this repo already has); the GENCODE
annotation parse itself is a local, already-cached file (no network),
so this just reads the committed demo cache, computes and embeds the
subset, and writes it back.
"""

from __future__ import annotations

from mudata import read_h5mu, write_h5mu

from multiome_agent.agent.demo_fixed_core_cache import DEMO_CACHE_PATH, _embed_gencode_gene_coords_subset
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def backfill() -> None:
    demo = read_h5mu(DEMO_CACHE_PATH)
    _embed_gencode_gene_coords_subset(demo)
    demo.update()
    write_h5mu(DEMO_CACHE_PATH, demo)
    logger.info("Backfilled %s with gencode_gene_coords_subset", DEMO_CACHE_PATH)


if __name__ == "__main__":
    backfill()
