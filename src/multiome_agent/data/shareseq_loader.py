"""Loader for the owner's shareseq-multi-cell-lines multi-condition multiome
data (see CLAUDE.md's "Data" section -- described only generically here and
everywhere else in this repo, never by file path, project, or lab name).

Three files, paths resolved from `config/local_paths.yaml` (gitignored,
never committed): an all-genes RNA AnnData, an HVG-only RNA AnnData (same
cells, precomputed highly-variable-gene subset), and an ATAC peak-count
AnnData. All three share the same ~5,800 cell barcodes in the same order
(same-cell multiome) -- verified programmatically below, not assumed.

Unlike the tenx-cell-ranger loader, these files need no download/subsample/cache
pipeline -- they're already single, reasonably-sized, pre-QC'd files. This
loader only aligns and combines them.
"""

from __future__ import annotations

import anndata as ad
from mudata import MuData

from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD
from multiome_agent.data.loader import _parse_peak_coords


def _require_path(path: str | None, config_key: str) -> str:
    if not path:
        raise RuntimeError(
            f"{config_key} is not set. Add it to the local, gitignored "
            "config/local_paths.yaml (see CLAUDE.md's Data section)."
        )
    return path


def load_shareseq_multiome() -> MuData:
    """Load and align the three shareseq files into one MuData.

    Returns a MuData with three aligned views over the same cells:
    ``"rna"`` (all genes), ``"rna_hvg"`` (the precomputed HVG subset -- used
    for embeddings/clustering per CLAUDE.md's guidance, so the existing HVG
    selection is reused rather than recomputed), and ``"atac"`` (raw peak
    counts). A third "modality" alongside rna/atac is unusual, but MuData
    supports it fine as long as obs align, and it avoids forcing two RNA
    views that share genes (not true here -- HVG is a *subset* of the
    all-genes var space) into a single AnnData with two layers.
    """
    rna_path = _require_path(SHARESEQ_RNA_H5AD, "SHARESEQ_RNA_H5AD")
    atac_path = _require_path(SHARESEQ_ATAC_H5AD, "SHARESEQ_ATAC_H5AD")
    hvg_path = _require_path(SHARESEQ_RNA_HVG_H5AD, "SHARESEQ_RNA_HVG_H5AD")

    rna = ad.read_h5ad(rna_path)
    rna_hvg = ad.read_h5ad(hvg_path)
    atac = ad.read_h5ad(atac_path)

    # `.var` ships with zero columns (no nearest-gene annotation) -- but
    # coordinates alone (chrom/start/end, parsed from var_names, same
    # "chr1:794759-795260" format as the tenx-cell-ranger data) are all that's
    # needed for genome-sequence-based work (motif deviations); nearest-gene
    # annotation is a separate, still-unbuilt thing (needed for gene
    # activity, not for this). Reuses the tenx-cell-ranger loader's parser rather than
    # duplicating it.
    atac.var = atac.var.join(_parse_peak_coords(atac.var_names))

    if not (list(rna.obs_names) == list(rna_hvg.obs_names) == list(atac.obs_names)):
        raise ValueError(
            "Shareseq RNA/RNA-HVG/ATAC barcodes are not identical and in the "
            "same order -- loader assumes same-cell multiome alignment "
            "(previously verified during inspection; re-verify before "
            "relying on it, don't assume it still holds if source files change)."
        )

    return MuData({"rna": rna, "rna_hvg": rna_hvg, "atac": atac})
