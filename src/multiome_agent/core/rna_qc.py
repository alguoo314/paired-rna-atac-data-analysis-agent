"""Fixed-core RNA QC: genes/UMIs per cell, % mito, doublet scores.

Operates on the RNA AnnData produced by ``multiome_agent.data.loader``
(raw counts, no existing QC columns) and writes results into ``.obs``.
"""

from __future__ import annotations

import scanpy as sc
from anndata import AnnData

from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def run_rna_qc(rna: AnnData) -> AnnData:
    """Compute RNA QC metrics in place and return the same object.

    Adds to ``.obs``: ``n_genes_by_counts``, ``total_counts``,
    ``pct_counts_mt``, ``doublet_score``, ``predicted_doublet``.
    """
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    n_mt = int(rna.var["mt"].sum())
    logger.info("Flagged %d mitochondrial genes (var_names starting with 'MT-')", n_mt)

    sc.pp.calculate_qc_metrics(
        rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # Scrublet needs raw counts and its own internal preprocessing; run on a
    # copy's `.X` so it doesn't fight with any later normalization the caller
    # does on `rna` itself. `copy=False` still mutates `rna` in place here
    # since we pass `rna` directly (no separate working copy needed at this
    # stage -- QC runs before any normalization/clustering).
    sc.pp.scrublet(rna)

    logger.info(
        "RNA QC: median genes/cell=%.0f, median UMIs/cell=%.0f, "
        "median pct_mt=%.2f, predicted doublets=%d/%d",
        rna.obs["n_genes_by_counts"].median(),
        rna.obs["total_counts"].median(),
        rna.obs["pct_counts_mt"].median(),
        int(rna.obs["predicted_doublet"].sum()),
        rna.n_obs,
    )
    return rna
