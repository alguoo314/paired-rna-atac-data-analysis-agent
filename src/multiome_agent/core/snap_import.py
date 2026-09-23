"""Shared snapatac2 fragments import, reused by ATAC QC (TSS enrichment) and
gene activity scoring so the ~2GB fragments file only gets sorted/parsed once
per pipeline run instead of twice.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import snapatac2 as snap

from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def import_snap_fragments(barcodes: Iterable[str], fragments_file: Path):
    """Import fragments for exactly our subsampled cells via snapatac2.

    Returns a snapatac2 AnnData indexed by cell barcode, with per-fragment
    insertions in `.obsm["insertion"]` (used downstream by both `tsse` and
    `make_gene_matrix`).
    """
    barcodes = list(barcodes)
    logger.info(
        "snapatac2.pp.import_fragments: %d whitelisted barcodes (sorts the "
        "whole fragments file regardless of whitelist size -- a few minutes)",
        len(barcodes),
    )
    data = snap.pp.import_fragments(
        fragments_file,
        chrom_sizes=snap.genome.hg38,
        whitelist=barcodes,
        min_num_fragments=0,
        sorted_by_barcode=False,
    )
    logger.info("snapatac2 import complete: %d cells", data.n_obs)
    return data
