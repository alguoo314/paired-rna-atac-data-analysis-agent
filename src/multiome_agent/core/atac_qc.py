"""Fixed-core ATAC QC: fragments/cell, FRiP, TSS enrichment, nucleosome signal.

Fragments-per-cell and FRiP come straight from Cell Ranger ARC's per-barcode
metrics (already in `.obs` via the loader). TSS enrichment and nucleosome
signal both need the real fragments file, via two different validated
wrappers from CLAUDE.md's approved stack:

- TSS enrichment: `snapatac2.metrics.tsse`, on a snapatac2 AnnData already
  imported (restricted to our subsampled barcodes) via
  `multiome_agent.core.snap_import.import_snap_fragments` -- shared with
  gene-activity scoring so the ~2GB fragments file is sorted/parsed only
  once per pipeline run, not once per consumer.
- Nucleosome signal: `muon.atac.tl.nucleosome_signal`, which samples
  fragments directly from the tabix-indexed file for just our cells (fast --
  seconds, not minutes) and returns the standard ArchR/Signac-style ratio of
  mono-nucleosomal (147-294bp) to nucleosome-free (<147bp) fragments.
  snapatac2 doesn't expose an equivalent *per-cell* metric (its
  `frag_size_distr` is dataset-level only), so muon is the right tool here.
"""

from __future__ import annotations

from pathlib import Path

import muon as mu
import snapatac2 as snap
from anndata import AnnData

from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def _compute_frip(atac: AnnData) -> None:
    atac.obs["frip"] = atac.obs["atac_peak_region_fragments"] / atac.obs["atac_fragments"]


def _compute_tss_enrichment(atac: AnnData, snap_data) -> None:
    snap.metrics.tsse(snap_data, snap.genome.hg38)
    tsse = snap_data.obs["tsse"].reindex(atac.obs_names)
    if tsse.isna().any():
        missing = int(tsse.isna().sum())
        raise ValueError(f"{missing} cell barcodes got no TSS enrichment score from snapatac2")
    atac.obs["tss_enrichment"] = tsse.to_numpy()


def _compute_nucleosome_signal(atac: AnnData, fragments_file: Path) -> None:
    mu.atac.tl.locate_fragments(atac, str(fragments_file))
    mu.atac.tl.nucleosome_signal(atac)


def run_atac_qc(atac: AnnData, fragments_file: Path, snap_data=None) -> AnnData:
    """Compute ATAC QC metrics in place and return the same object.

    Adds to ``.obs``: ``frip``, ``tss_enrichment``, ``nucleosome_signal``
    (``atac_fragments`` already present from the loader).

    `snap_data`: an already-imported snapatac2 AnnData for these same cells
    (see `snap_import.import_snap_fragments`) -- pass it in if the caller
    already built one (e.g. for gene activity too), otherwise one is built
    here from `fragments_file`.
    """
    if snap_data is None:
        from multiome_agent.core.snap_import import import_snap_fragments

        snap_data = import_snap_fragments(atac.obs_names, fragments_file)

    _compute_frip(atac)
    _compute_tss_enrichment(atac, snap_data)
    _compute_nucleosome_signal(atac, fragments_file)

    logger.info(
        "ATAC QC: median FRiP=%.3f, median TSS enrichment=%.2f, median nucleosome signal=%.2f",
        atac.obs["frip"].median(),
        atac.obs["tss_enrichment"].median(),
        atac.obs["nucleosome_signal"].median(),
    )
    return atac
