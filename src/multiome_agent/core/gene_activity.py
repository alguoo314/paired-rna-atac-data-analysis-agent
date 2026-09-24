"""Fixed-core gene activity scores: ATAC accessibility near each gene as a
proxy for gene activity, via `snapatac2.pp.make_gene_matrix` (counts TN5
insertions in each gene's TSS + gene body, extended by an up/downstream
window -- the standard ArchR/Signac-style gene activity definition).

Needs a snapatac2 AnnData already imported from the fragments file (see
`multiome_agent.core.snap_import.import_snap_fragments`) -- shared with ATAC
QC's TSS enrichment step so the fragments file is only imported once.

`compute_gene_activity_from_peaks` below is a second, independent method for
datasets that have a peak-count matrix but no fragments file (see its
docstring) -- it needs no fragments import at all.
"""

from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp
import snapatac2 as snap
from mudata import MuData


def compute_gene_activity(mdata: MuData, snap_data) -> pd.DataFrame:
    """Compute a cells x genes gene-activity matrix and attach it to `mdata`.

    Stored as a raw scipy sparse matrix in
    `mdata.mod["atac"].obsm["gene_activity"]` plus the gene names in
    `mdata.mod["atac"].uns["gene_activity_genes"]` (mirroring how
    `chromvar_deviations`/`chromvar_motif_names` are split in
    motif_deviations.py) -- NOT a `pandas.DataFrame.sparse` (pandas' sparse
    *extension dtype*, not scipy sparse): anndata's h5mu writer has no
    registered method for that pandas dtype and raises `IORegistryError` on
    write, only surfacing when something actually caches the MuData to disk
    (not in step 4's in-memory-only tests, only once step 5's agent cache
    tried to persist it). A plain scipy sparse matrix in `.obsm` is natively
    h5-writable. It's a derived *view* of ATAC accessibility keyed by gene
    rather than peak, not an independently-measured modality, so it belongs
    alongside the ATAC AnnData rather than sitting next to "rna"/"atac" as a
    peer modality.
    """
    atac = mdata.mod["atac"]
    gene_mat = snap.pp.make_gene_matrix(snap_data, gene_anno=snap.genome.hg38, inplace=False)

    # snap_data's cell order need not match `atac`'s -- reindex explicitly.
    gene_mat = gene_mat[atac.obs_names].copy()
    atac.obsm["gene_activity"] = gene_mat.X
    atac.uns["gene_activity_genes"] = list(gene_mat.var_names)
    return pd.DataFrame.sparse.from_spmatrix(
        gene_mat.X, index=gene_mat.obs_names, columns=gene_mat.var_names
    )


def _load_protein_coding_gene_windows(upstream_bp: int = 2000) -> pd.DataFrame:
    """Gene coordinate windows for peak-summation-based gene activity (see
    `compute_gene_activity_from_peaks`), parsed directly from the same
    cached GENCODE annotation `compute_gene_activity` already uses
    (`snap.genome.hg38.annotation`, a local GFF3 file snapatac2 downloads
    once and reuses). Reusable as-is for ANY hg38/GRCh38 dataset -- human
    gene coordinates are a property of the genome build, not of which cells
    were sequenced -- so this needs no separate download for a new dataset,
    only a genome-build check that the dataset's own peaks are really hg38
    (done by the caller, not here; see `shareseq_pipeline._sanity_check_peak_sequences`).

    Restricted to `gene_type=protein_coding` (~20K of GENCODE's ~62K "gene"
    rows) since those are what marker/TF lookups actually need, not
    pseudogenes/lncRNAs. A handful of gene symbols (7 among autosomal
    protein-coding genes, e.g. duplicated paralog annotations like
    HERC3/MATR3) appear as more than one GENCODE "gene" record; these are
    merged into a single window per symbol (min start, max end) rather than
    picked arbitrarily, checked directly against this exact file rather
    than assumed.
    """
    path = snap.genome.hg38.annotation
    rows = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene" or "gene_type=protein_coding" not in fields[8]:
                continue
            match = re.search(r"gene_name=([^;]+)", fields[8])
            if match:
                rows.append((fields[0], int(fields[3]), int(fields[4]), fields[6], match.group(1)))

    genes = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene_name"])
    genes = genes.groupby(["chrom", "gene_name", "strand"], as_index=False).agg(
        start=("start", "min"), end=("end", "max")
    )
    # Promoter-extend upstream of the TSS (start for "+", end for "-"),
    # leaving the downstream gene-body boundary untouched -- ArchR/Signac's
    # default gene-activity window.
    win_start = np.where(genes["strand"] == "+", genes["start"] - upstream_bp, genes["start"])
    win_end = np.where(genes["strand"] == "+", genes["end"], genes["end"] + upstream_bp)
    genes["win_start"] = np.clip(win_start, 0, None)
    genes["win_end"] = win_end
    return genes[["chrom", "gene_name", "win_start", "win_end"]]


def _sum_peaks_into_gene_windows(
    atac_X_csc, peaks: pd.DataFrame, gene_windows: pd.DataFrame
) -> tuple[sp.csr_matrix, list[str]]:
    """Pure overlap-and-sum core of `compute_gene_activity_from_peaks`, kept
    separate from genome-annotation loading so it's unit-testable on small
    synthetic coordinate tables without touching the real GENCODE file.

    `peaks` (chrom/start/end) must be in the same row order as `atac_X_csc`'s
    columns. A peak counts toward a gene if its interval overlaps that
    gene's window at all (interval overlap, not point containment) -- one
    peak can contribute to more than one gene when nearby genes' windows
    overlap, matching ArchR/Signac's definition rather than picking a single
    "nearest gene" per peak.
    """
    peaks = peaks.reset_index(drop=True)
    gene_names: list[str] = []
    gene_cols = []
    for chrom, gsub in gene_windows.groupby("chrom"):
        psub = peaks[peaks["chrom"] == chrom]
        if psub.empty or gsub.empty:
            continue
        p_start, p_end = psub["start"].to_numpy(), psub["end"].to_numpy()
        p_positions = psub.index.to_numpy()
        g_wstart, g_wend = gsub["win_start"].to_numpy(), gsub["win_end"].to_numpy()
        # genes x peaks overlap matrix via broadcasting, then peaks x genes
        # for the matmul below.
        overlap = (p_start[None, :] < g_wend[:, None]) & (p_end[None, :] > g_wstart[:, None])
        if not overlap.any():
            continue
        overlap_sp = sp.csr_matrix(overlap.astype(np.float64))
        gene_cols.append(atac_X_csc[:, p_positions] @ overlap_sp.T)
        gene_names.extend(gsub["gene_name"].tolist())

    if not gene_cols:
        return sp.csr_matrix((atac_X_csc.shape[0], 0)), []
    return sp.hstack(gene_cols).tocsr(), gene_names


def compute_gene_activity_from_peaks(mdata: MuData, upstream_bp: int = 2000) -> pd.DataFrame:
    """Alternative gene-activity method for ATAC data that has a peak-count
    matrix but no fragments file (the shareseq-multi-cell-lines dataset --
    `compute_gene_activity`'s `snap.pp.make_gene_matrix` needs a snapatac2-
    imported fragments object, which doesn't exist for that data; see
    PROGRESS.md's step 19 entry on this).

    This is the standard fallback ArchR/Signac use when only a peak matrix
    is available: sum each cell's peak counts over whichever peaks overlap
    a gene's window (gene body + an `upstream_bp` promoter extension). Uses
    the SAME cached GENCODE annotation `compute_gene_activity` uses --
    reusable here because gene coordinates depend only on genome build, not
    on which cells were sequenced, and callers are expected to have already
    confirmed this dataset's peaks are really hg38 (see
    `_load_protein_coding_gene_windows`'s docstring).

    Stores results the same way as `compute_gene_activity`
    (`atac.obsm["gene_activity"]` / `atac.uns["gene_activity_genes"]`) so
    downstream code doesn't need to know which method produced them.
    """
    atac = mdata.mod["atac"]
    gene_windows = _load_protein_coding_gene_windows(upstream_bp)
    gene_activity, gene_names = _sum_peaks_into_gene_windows(
        atac.X.tocsc(), atac.var[["chrom", "start", "end"]], gene_windows
    )
    atac.obsm["gene_activity"] = gene_activity
    atac.uns["gene_activity_genes"] = gene_names
    return pd.DataFrame.sparse.from_spmatrix(
        gene_activity, index=atac.obs_names, columns=gene_names
    )
