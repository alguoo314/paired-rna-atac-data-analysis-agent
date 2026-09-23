"""Corrupt copies of the clean PBMC dataset for the fault-injection benchmark.

CLAUDE.md's fault-injection spec lists five fault types: cell-line/cell-type
label swaps, shuffled RNA-ATAC barcode pairing, heavy downsampling of one
modality, mixed samples, and injected doublets. This is a single-sample PBMC
dev dataset with no cell-line/condition labels, so "cell-line label swaps"
and "mixed samples" don't apply here -- both need the multi-cell-line
showcase data (a later stretch goal, not built). Three fault types are
implemented, each with severity levels:

1. Shuffled RNA-ATAC pairing (`shuffle_rna_atac_pairing`) -- permutes which
   ATAC profile sits under which RNA barcode. Per-cell computations
   (each modality's own QC, gene activity, motif deviations) are UNCHANGED
   by this fault, since they don't depend on cross-modal pairing -- only
   cross-modal quantities (cluster agreement, TF-motif correlation) shift.
   This means a fixed-core-processed clean result can be reused directly
   (just permute+relabel the ATAC AnnData's rows and recompute the two
   cross-modal functions), with no expensive snapatac2 rerun needed.
2. Heavy ATAC downsampling (`downsample_atac_fragments`) -- genuinely drops
   a fraction of fragment lines before snapatac2 import, so the fault
   propagates into every fragments-derived metric (TSS enrichment,
   nucleosome signal, gene activity), not just peak counts. This DOES need
   a fresh snapatac2 import on the downsampled fragments file.
3. Injected doublets (`inject_doublets`) -- synthesizes new cells by summing
   raw RNA counts and raw ATAC peak counts of random real-cell pairs (a
   real doublet is one droplet contributing to both readouts). Synthetic
   cells have no real fragments in the source file, so fragments-derived
   ATAC QC (TSS enrichment, nucleosome signal) is left as NaN for them
   rather than faked -- see `inject_doublets`'s docstring.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pysam
import scipy.sparse as sp
from anndata import AnnData
from mudata import MuData


@dataclass
class FaultMetadata:
    fault_type: str
    severity: float
    description: str
    ground_truth: dict = field(default_factory=dict)


def clean_control(mdata: MuData) -> tuple[MuData, FaultMetadata]:
    """Unmodified passthrough, for false-alarm measurement."""
    return mdata.copy(), FaultMetadata(
        fault_type="clean_control", severity=0.0,
        description="Unmodified data (no fault injected).",
    )


def shuffle_rna_atac_pairing(mdata: MuData, fraction: float, seed: int = 0) -> tuple[MuData, FaultMetadata]:
    """Permute which ATAC profile is paired with which RNA barcode.

    `fraction` of cells (by index, not by any biological property) have
    their ATAC profile swapped with another shuffled cell's, via a
    fixed-point-free permutation restricted to the selected subset (so
    "shuffled" cells are guaranteed to actually change pairing, not
    coincidentally map to themselves). The remaining `1 - fraction` keep
    their true pairing. RNA is untouched; only the ATAC AnnData's row order
    is permuted and then relabeled back to the original obs_names, so each
    barcode's per-cell computed values (QC, gene activity, motif deviations)
    travel WITH their true source cell -- correct for the "reuse a clean
    fixed-core result" optimization used by the eval scoring harness.
    """
    rng = np.random.default_rng(seed)
    atac = mdata.mod["atac"].copy()
    n = atac.n_obs
    obs_names = list(atac.obs_names)

    n_shuffle = int(round(fraction * n))
    shuffle_idx = rng.choice(n, size=n_shuffle, replace=False)
    perm = np.arange(n)
    if n_shuffle > 1:
        sub = shuffle_idx.copy()
        # Fixed-point-free permutation of `sub`: rotate by a random non-zero
        # offset. A rotation can't fix any element as long as the offset
        # isn't a multiple of len(sub), guaranteed here since 0 < offset < len(sub).
        offset = rng.integers(1, len(sub))
        perm[shuffle_idx] = np.roll(sub, offset)
    elif n_shuffle == 1:
        # A single selected index has no valid derangement partner; drop it
        # rather than silently leave it unshuffled (fraction would then not
        # match the actual shuffled count).
        n_shuffle = 0

    atac_shuffled = atac[perm].copy()
    atac_shuffled.obs_names = obs_names

    true_source = {obs_names[i]: obs_names[perm[i]] for i in range(n) if perm[i] != i}

    mdata_new = MuData({"rna": mdata.mod["rna"].copy(), "atac": atac_shuffled})
    mdata_new.uns.update(mdata.uns)

    return mdata_new, FaultMetadata(
        fault_type="shuffled_rna_atac_pairing",
        severity=fraction,
        description=f"{len(true_source)}/{n} cells have their ATAC profile swapped with another cell's.",
        ground_truth={"true_atac_source_by_rna_barcode": true_source},
    )


def _stream_downsample_fragments(
    fragments_file: Path, barcodes: set[str], keep_frac: float, seed: int, out_path: Path
) -> int:
    """Write a new PLAIN-TEXT fragments file (not gzipped -- `pysam.tabix_index`
    below does the bgzip compression itself) containing only lines whose
    barcode is in `barcodes`, each independently kept with probability
    `keep_frac`. Preserves the input's genomic sort order (we only drop
    lines, never reorder), which `tabix_index` requires. Returns the number
    of lines written."""
    rng = np.random.default_rng(seed)
    n_written = 0
    with gzip.open(fragments_file, "rt") as fin, open(out_path, "w") as fout:
        for line in fin:
            barcode = line.split("\t", 4)[3]
            if barcode in barcodes and rng.random() < keep_frac:
                fout.write(line)
                n_written += 1
    return n_written


def downsample_atac_fragments(
    fragments_file: Path, barcodes: list[str], keep_frac: float, out_dir: Path, seed: int = 0
) -> tuple[Path, FaultMetadata]:
    """Write a downsampled copy of `fragments_file`, restricted to `barcodes`
    and keeping `keep_frac` of their fragment lines, then bgzip-compress and
    tabix-index it (via `pysam.tabix_index`) so it's usable both by
    `snap.pp.import_fragments` (which sorts internally so doesn't strictly
    need the index, but accepts it fine) AND by `muon.atac.tl.locate_fragments`
    (which uses `pysam.TabixFile` and silently fails -- caught by muon's own
    broad `except Exception: print(...)`, no exception raised, leaving
    whatever fragments path was previously set in `.uns` -- if the file
    ISN'T tabix-indexed; a real bug caught by inspecting `nucleosome_signal`
    output for a fault that should have shifted it but silently didn't).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    plain_path = out_dir / f"downsampled_keep{keep_frac:.2f}_seed{seed}.tsv"
    n_written = _stream_downsample_fragments(fragments_file, set(barcodes), keep_frac, seed, plain_path)
    out_path = Path(pysam.tabix_index(str(plain_path), preset="bed", force=True))

    return out_path, FaultMetadata(
        fault_type="atac_downsampling",
        severity=1.0 - keep_frac,
        description=f"ATAC fragments downsampled to {keep_frac:.0%} (kept {n_written} lines for {len(barcodes)} barcodes).",
        ground_truth={"keep_frac": keep_frac, "n_fragments_kept": n_written},
    )


def inject_doublets(mdata: MuData, doublet_rate: float, seed: int = 0) -> tuple[MuData, FaultMetadata]:
    """Append synthetic doublet cells (summed raw RNA + ATAC counts of
    random real-cell pairs) to `mdata`. `doublet_rate` is the fraction of
    the FINAL cell count that is synthetic (e.g. 0.1 on 500 real cells adds
    ~56 synthetic cells so synthetic/(real+synthetic) ~= 0.1).

    Expects `mdata` to carry RAW counts (as `load_pbmc_multiome` produces,
    before `run_fixed_core`'s normalization) -- summing log-normalized
    values would not correspond to a real doublet's physical count mixing.

    Synthetic cells have no entry in the source fragments file (they're not
    real droplets), so fragments-derived ATAC QC (`tss_enrichment`,
    `nucleosome_signal`, and `atac_fragments`-based `frip`) cannot be
    computed for them -- callers must handle/exclude synthetic cells for
    those specific fields rather than expect real values.
    """
    rng = np.random.default_rng(seed)
    rna = mdata.mod["rna"].copy()
    atac = mdata.mod["atac"].copy()
    n_real = rna.n_obs

    n_doublets = int(round(doublet_rate * n_real / (1 - doublet_rate))) if doublet_rate < 1 else n_real
    pair_a = rng.integers(0, n_real, size=n_doublets)
    pair_b = rng.integers(0, n_real, size=n_doublets)
    # A cell paired with itself isn't a doublet (no new information); redraw.
    same = pair_a == pair_b
    while same.any():
        pair_b[same] = rng.integers(0, n_real, size=same.sum())
        same = pair_a == pair_b

    doublet_names = [f"synthetic_doublet_{i:05d}" for i in range(n_doublets)]

    def _append_summed_pairs(adata: AnnData) -> AnnData:
        summed = adata.X[pair_a] + adata.X[pair_b]
        new_X = sp.vstack([adata.X, summed]).tocsr() if sp.issparse(adata.X) else np.vstack([adata.X, summed])
        new_obs = adata.obs.iloc[pair_a].copy()
        new_obs.index = doublet_names
        return ad.AnnData(X=new_X, obs=pd.concat([adata.obs, new_obs], axis=0), var=adata.var.copy())

    rna_new = _append_summed_pairs(rna)
    atac_new = _append_summed_pairs(atac)
    for col in ("atac_fragments", "atac_TSS_fragments", "atac_peak_region_fragments"):
        if col in atac_new.obs.columns:
            atac_new.obs.loc[doublet_names, col] = np.nan
    for col in ("gex_umis_count", "gex_genes_count"):
        if col in rna_new.obs.columns:
            rna_new.obs.loc[doublet_names, col] = np.nan

    mdata_new = MuData({"rna": rna_new, "atac": atac_new})
    mdata_new.uns.update(mdata.uns)

    return mdata_new, FaultMetadata(
        fault_type="injected_doublets",
        severity=doublet_rate,
        description=f"{n_doublets} synthetic doublets added to {n_real} real cells ({n_doublets / (n_real + n_doublets):.1%} of final total).",
        ground_truth={"doublet_obs_names": doublet_names, "source_pairs": list(zip(pair_a.tolist(), pair_b.tolist()))},
    )
