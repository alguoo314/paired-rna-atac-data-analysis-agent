"""One-time repair for the already-committed demo cache
(`demo_cache/agent_demo_core.h5mu`), backfilling the fields
`_trim_to_tool_readable_fields` used to strip before `peak_to_gene_links`/
`regulon_inference` existed (see `agent/demo_fixed_core_cache.py`'s module
docstring for the full story): real peak coordinates (`chrom`/`start`/
`end`), the real per-cell peak-count matrix, and the peak x motif match
matrix.

Deliberately NOT a rebuild from raw files (that pays `build_demo_cache`'s
own ~17-20 min fragments-reimport cost again, for data this repo already
has computed). Instead, reindexes directly from the already-built, already
`motif_match`-persisted real agent fixed-core cache
(`data/processed/agent_fixed_core_n11909_seed0.h5mu`) -- verified before
writing this script, not assumed, that this is safe:
- the demo's 150 cell barcodes are a real subset of the full cache's
  11,909 (both are seeded draws over the exact same underlying PBMC
  dataset, just at different `n_cells`, with no guarantee of a NESTED
  draw -- confirmed by actually comparing `obs_names` sets, not assumed
  from the shared seed).
- the demo's own (independently re-filtered-for-150-cells, hence smaller:
  61,474 vs. 107,385) ATAC peak set is a real name-subset of the full
  cache's peak set -- confirmed the same way.

Run once; the output already matches `_trim_to_tool_readable_fields`'s
current (fixed) behavior, so a future full `build_demo_cache()` rebuild
from raw files would reproduce an equivalent object on its own.
"""

from __future__ import annotations

from mudata import read_h5mu, write_h5mu

from multiome_agent.agent.demo_fixed_core_cache import DEMO_CACHE_PATH
from multiome_agent.agent.fixed_core_cache import CACHE_PATH as FULL_CACHE_PATH
from multiome_agent.core.motif_deviations import motif_name_tokens
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def backfill() -> None:
    demo = read_h5mu(DEMO_CACHE_PATH)
    full = read_h5mu(FULL_CACHE_PATH)

    demo_atac, full_atac = demo.mod["atac"], full.mod["atac"]
    demo_obs, demo_peaks = list(demo_atac.obs_names), list(demo_atac.var_names)

    assert set(demo_obs).issubset(set(full_atac.obs_names)), "demo cells must be a subset of the full cache's cells"
    assert set(demo_peaks).issubset(set(full_atac.var_names)), "demo peaks must be a subset of the full cache's peaks"

    full_reindexed = full_atac[demo_obs, demo_peaks].copy()

    demo_atac.var["chrom"] = full_reindexed.var["chrom"].to_numpy()
    demo_atac.var["start"] = full_reindexed.var["start"].to_numpy()
    demo_atac.var["end"] = full_reindexed.var["end"].to_numpy()
    demo_atac.layers["counts"] = full_reindexed.layers["counts"]

    # Same filter `_trim_to_tool_readable_fields` now applies on a fresh build:
    # `regulon_inference` rejects any TF not in `rna.var_names` at its very first
    # gate regardless of whether its motif column exists, so dropping columns for
    # TFs absent from this RNA panel changes zero real behavior while cutting the
    # matrix from 879 columns down to the ones that matter -- the dominant
    # contributor to file size at this matrix's real density (sparse CSR costs
    # MORE than dense here, verified by measuring the actual h5 layout).
    rna_genes = {g.upper() for g in demo.mod["rna"].var_names}
    motif_names = list(full_atac.uns["motif_match_names"])
    keep_idx = [i for i, name in enumerate(motif_names) if any(t.upper() in rna_genes for t in motif_name_tokens(name))]
    demo_atac.varm["motif_match"] = full_reindexed.varm["motif_match"][:, keep_idx]
    demo_atac.uns["motif_match_names"] = [motif_names[i] for i in keep_idx]

    demo.update()
    write_h5mu(DEMO_CACHE_PATH, demo)
    logger.info(
        "Backfilled %s: chrom/start/end + motif_match (%s) + counts layer for %d cells x %d peaks",
        DEMO_CACHE_PATH, demo_atac.varm["motif_match"].shape, len(demo_obs), len(demo_peaks),
    )


if __name__ == "__main__":
    backfill()
