"""Code-level scoring: does a fault actually shift the fixed-core numbers?

This is deliberately independent of any LLM agent -- it answers "is this
fault detectable in the computed tables at all," which is a prerequisite
for an agent ever being able to notice it. Each `score_*` function compares
a faulted result against a clean baseline and returns a small delta dict
(no pass/fail judgment baked in -- that's the eval harness's job later).
"""

from __future__ import annotations

from pathlib import Path

import muon as mu
import numpy as np
import pandas as pd
import snapatac2 as snap
from anndata import AnnData
from mudata import MuData

from multiome_agent.core.clustering import ATAC_CLUSTER_KEY, RNA_CLUSTER_KEY, rna_atac_cluster_agreement
from multiome_agent.core.motif_deviations import tf_expression_motif_correlation
from multiome_agent.core.snap_import import import_snap_fragments
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def score_shuffled_pairing(clean_mdata: MuData, shuffled_mdata: MuData) -> dict:
    """Cross-modal signals only -- both cheap to recompute, no snapatac2 rerun.

    Per-cell signals (QC, gene activity, motif deviations) are identical by
    construction (shuffling only changes which ATAC row sits under which
    RNA barcode label) and are not re-scored here.
    """
    clean_agreement = clean_mdata.uns["rna_atac_cluster_agreement"]
    shuffled_agreement = rna_atac_cluster_agreement(shuffled_mdata)

    clean_tf = clean_mdata.uns["tf_motif_sanity_check"]
    shuffled_tf = tf_expression_motif_correlation(shuffled_mdata, gene=clean_tf["gene"], motif_name_contains=clean_tf["gene"])

    return {
        "ari_clean": clean_agreement["ari"],
        "ari_faulted": shuffled_agreement["ari"],
        "ari_delta": shuffled_agreement["ari"] - clean_agreement["ari"],
        "tf_motif_rho_clean": clean_tf["spearman_rho"],
        "tf_motif_rho_faulted": shuffled_tf["spearman_rho"],
        "tf_motif_rho_delta": shuffled_tf["spearman_rho"] - clean_tf["spearman_rho"],
    }


def _downsampled_atac_metrics(atac_barcodes: list[str], fragments_path: Path, clean_atac: AnnData) -> pd.DataFrame:
    """n_fragment / TSS enrichment / nucleosome signal directly from a
    (possibly downsampled) fragments file -- NOT `atac_qc.run_atac_qc`,
    whose `frip` reuses Cell Ranger's original (un-downsampled)
    `atac_peak_region_fragments`/`atac_fragments` metrics from
    `per_barcode_metrics.csv`, which wouldn't reflect this fault. Fragment
    count, TSS enrichment, and nucleosome signal ARE recomputed directly
    from the (possibly downsampled) fragments file, so they genuinely
    reflect it.
    """
    snap_data = import_snap_fragments(atac_barcodes, fragments_path)
    snap.metrics.tsse(snap_data, snap.genome.hg38)

    n_fragment = snap_data.obs["n_fragment"].reindex(atac_barcodes)

    probe = clean_atac[atac_barcodes].copy()
    mu.atac.tl.locate_fragments(probe, str(fragments_path))
    # muon's `nucleosome_signal` does NOT iterate per barcode -- it does one
    # genome-wide `fragments.fetch()` and reads a fixed `n` fragments total
    # (default `adata.n_obs * 1e4`), bucketing each by barcode as it goes; its
    # `except KeyError: pass` only swallows unknown-barcode lines, not the
    # `StopIteration` `next(fr)` raises once the file is exhausted. A
    # downsampled (or barcode-restricted) file can easily have fewer total
    # lines than that fixed target -- e.g. keep_frac=0.2 keeps ~20% of lines,
    # so the default `n` massively overshoots -- so it crashes the whole
    # batch rather than just stopping early. Cap `n` at the real total
    # fragment count in this exact file (already known from the snapatac2
    # import above) so it never tries to read past the end.
    total_fragments = int(n_fragment.sum())
    default_n = len(atac_barcodes) * 10_000
    mu.atac.tl.nucleosome_signal(probe, n=min(default_n, total_fragments))

    return pd.DataFrame({
        "n_fragment": n_fragment,
        "tss_enrichment": snap_data.obs["tsse"].reindex(atac_barcodes),
        "nucleosome_signal": probe.obs["nucleosome_signal"].reindex(atac_barcodes),
    })


def score_downsampling(
    clean_atac: AnnData, atac_barcodes: list[str], clean_fragments_path: Path, downsampled_fragments_path: Path, keep_frac: float
) -> dict:
    """Compare fragments-derived ATAC metrics between the full fragments
    file (restricted to `atac_barcodes`) and a downsampled copy.
    """
    clean_metrics = _downsampled_atac_metrics(atac_barcodes, clean_fragments_path, clean_atac)
    faulted_metrics = _downsampled_atac_metrics(atac_barcodes, downsampled_fragments_path, clean_atac)

    result = {"keep_frac": keep_frac}
    for col in clean_metrics.columns:
        c, f = clean_metrics[col].median(), faulted_metrics[col].median()
        result[f"{col}_median_clean"] = c
        result[f"{col}_median_faulted"] = f
        result[f"{col}_median_pct_change"] = (f - c) / c * 100 if c else float("nan")
    return result


def score_doublets(rna_obs_faulted: pd.DataFrame, doublet_obs_names: list[str]) -> dict:
    """Is the doublet-detection rate elevated among the synthetic cells vs.
    the real ones, after rerunning RNA QC (incl. scrublet) on the augmented
    dataset? This is the code-level "is this fault detectable" signal --
    independent of whether an LLM agent later notices it.
    """
    is_synthetic = rna_obs_faulted.index.isin(doublet_obs_names)
    real_rate = rna_obs_faulted.loc[~is_synthetic, "predicted_doublet"].mean()
    synthetic_rate = rna_obs_faulted.loc[is_synthetic, "predicted_doublet"].mean()
    real_score = rna_obs_faulted.loc[~is_synthetic, "doublet_score"].median()
    synthetic_score = rna_obs_faulted.loc[is_synthetic, "doublet_score"].median()

    return {
        "n_real": int((~is_synthetic).sum()),
        "n_synthetic": int(is_synthetic.sum()),
        "predicted_doublet_rate_real": float(real_rate),
        "predicted_doublet_rate_synthetic": float(synthetic_rate),
        "doublet_score_median_real": float(real_score),
        "doublet_score_median_synthetic": float(synthetic_score),
    }
