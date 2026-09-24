"""Load, subsample, and cache the 10x PBMC 10k Multiome dev dataset.

Inspection notes from the raw Cell Ranger ARC output (recorded here so later
readers don't have to re-derive them):

- ``filtered_feature_bc_matrix.h5`` is a single combined matrix: 11909 cells x
  144978 features. ``var['feature_types']`` is either ``"Gene Expression"``
  (36601 genes) or ``"Peaks"`` (108377 peaks) -- one file, two modalities, same
  cell barcodes (this is same-cell multiome, so RNA and ATAC are paired by
  construction; no separate barcode-matching step is needed). Peak
  ``var_names`` look like ``"chr1:10109-10357"``. Counts are raw integers
  (``float32``-typed but integer-valued), no ``.layers``. Single sample/genome
  (GRCh38) -- there's no cell-line or condition column for this dataset, only
  the shareseq-multi-cell-lines data (see CLAUDE.md) has those.
- ``per_barcode_metrics.csv`` has one row per *tested* barcode (732271 rows,
  most of which are empty droplets), with an ``is_cell`` flag; exactly 11909
  rows have ``is_cell == 1``, matching the filtered matrix's cell count
  exactly. Its ``barcode`` column matches ``adata.obs_names`` format exactly
  (e.g. ``AAACAGCCAAGGAATC-1``).
- ``atac_peak_annotation.tsv`` has exactly one row per peak (108377, no
  duplicates) keyed by ``peak`` in ``"chr1_10109_10357"`` format (underscores,
  not the colon/dash format used in ``var_names`` -- the two need reformatting
  to join). ``gene``/``distance``/``peak_type`` are ``";"``-joined when a peak
  sits near more than one gene (28% of peaks); we keep the raw joined string
  plus a convenience ``primary_gene`` (first token).

``atac_fragments.tsv.gz``/``.tbi`` and the hg38 reference genome are untouched
here -- they're for steps 3-4 (ATAC QC via fragments, motif deviations),
not this loader.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from mudata import MuData, read_h5mu, write_h5mu

from multiome_agent.config import REPO_ROOT, TENX_MATRIX_H5
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

RAW_DIR = REPO_ROOT / "data" / "raw" / "pbmc10k_multiome"
CACHE_DIR = REPO_ROOT / "data" / "processed"

def _resolve_matrix_h5(override: str | None, raw_dir: Path) -> Path:
    """`config/local_paths.yaml`'s `tenx_matrix_h5` key overrides the
    default -- lets a user point the combined-single-file loading strategy
    at their own 10x Cell Ranger ARC .h5 file without editing code (real
    limitation flagged in the README, now fixed: this used to be an
    unconditionally hardcoded path). Factored out so the branch is directly
    unit-testable without needing to reload this module after monkeypatching
    config at import time.
    """
    return Path(override) if override else raw_dir / "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"


MATRIX_H5 = _resolve_matrix_h5(TENX_MATRIX_H5, RAW_DIR)
PER_BARCODE_METRICS_CSV = RAW_DIR / "pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv"
PEAK_ANNOTATION_TSV = RAW_DIR / "pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv"

# Per-cell columns from Cell Ranger ARC worth carrying into .obs for later QC
# steps (step 3). Not exhaustive -- just the ones the fixed-core QC will need.
BARCODE_METRIC_COLUMNS = [
    "gex_umis_count",
    "gex_genes_count",
    "atac_fragments",
    "atac_TSS_fragments",
    "atac_peak_region_fragments",
]


def _parse_peak_coords(var_names: pd.Index) -> pd.DataFrame:
    """Parse "chr1:10109-10357"-style peak names into chrom/start/end columns."""
    names = var_names.to_series()
    chrom_coords = names.str.split(":", n=1, expand=True)
    start_end = chrom_coords[1].str.split("-", n=1, expand=True)
    return pd.DataFrame(
        {
            "chrom": chrom_coords[0],
            "start": start_end[0].astype(int),
            "end": start_end[1].astype(int),
        },
        index=var_names,
    )


def _split_rna_atac(adata: AnnData) -> tuple[AnnData, AnnData]:
    is_rna = adata.var["feature_types"] == "Gene Expression"
    is_atac = adata.var["feature_types"] == "Peaks"
    if not is_rna.any() or not is_atac.any():
        raise ValueError(
            f"Expected both 'Gene Expression' and 'Peaks' feature types, "
            f"got: {adata.var['feature_types'].unique().tolist()}"
        )
    rna = adata[:, is_rna].copy()
    atac = adata[:, is_atac].copy()
    return rna, atac


def _attach_barcode_metrics(obs: pd.DataFrame, metrics_csv=PER_BARCODE_METRICS_CSV) -> pd.DataFrame:
    metrics = pd.read_csv(metrics_csv)
    metrics = metrics[metrics["is_cell"] == 1].set_index("barcode")
    missing = obs.index.difference(metrics.index)
    if len(missing) > 0:
        raise ValueError(f"{len(missing)} cell barcodes missing from {metrics_csv.name}")
    return obs.join(metrics.loc[obs.index, BARCODE_METRIC_COLUMNS])


def _attach_peak_annotation(var: pd.DataFrame, annotation_tsv=PEAK_ANNOTATION_TSV) -> pd.DataFrame:
    ann = pd.read_csv(annotation_tsv, sep="\t")
    # annotation uses "chr1_10109_10357"; var_names use "chr1:10109-10357" -- reformat to join.
    ann_key = ann["peak"].str.replace("_", ":", n=1).str.replace("_", "-", n=1)
    ann = ann.set_index(ann_key).rename(
        columns={"gene": "nearest_genes", "peak_type": "peak_types"}
    )
    var = var.join(ann[["nearest_genes", "distance", "peak_types"]])
    var["primary_gene"] = var["nearest_genes"].str.split(";").str[0]
    return var


def _load_raw() -> MuData:
    """Load the raw Cell Ranger ARC h5 + side tables into a paired MuData."""
    logger.info("Loading raw matrix from %s", MATRIX_H5)
    adata = sc.read_10x_h5(MATRIX_H5, gex_only=False)
    adata.var_names_make_unique()

    rna, atac = _split_rna_atac(adata)

    peak_coords = _parse_peak_coords(atac.var_names)
    atac.var = atac.var.join(peak_coords)
    atac.var = _attach_peak_annotation(atac.var)

    rna.obs = _attach_barcode_metrics(rna.obs)
    atac.obs = _attach_barcode_metrics(atac.obs)

    mdata = MuData({"rna": rna, "atac": atac})
    logger.info(
        "Loaded MuData: %d cells, rna=%d genes, atac=%d peaks",
        mdata.n_obs, rna.n_vars, atac.n_vars,
    )
    return mdata


def subsample_mudata(mdata: MuData, n_cells: int = 3000, seed: int = 0) -> MuData:
    """Uniformly subsample cells (paired across modalities) for fast iteration."""
    n_cells = min(n_cells, mdata.n_obs)
    rng = np.random.default_rng(seed)
    keep = rng.choice(mdata.obs_names, size=n_cells, replace=False)
    sub = mdata[keep].copy()
    logger.info("Subsampled %d -> %d cells (seed=%d)", mdata.n_obs, sub.n_obs, seed)
    return sub


def _cache_path(n_cells: int, seed: int) -> Path:
    return CACHE_DIR / f"pbmc10k_multiome_n{n_cells}_seed{seed}.h5mu"


def load_pbmc_multiome(n_cells: int = 3000, seed: int = 0, force_reload: bool = False) -> MuData:
    """Load a subsampled, paired PBMC RNA+ATAC MuData, using a local cache.

    First call builds the subsample from the raw Cell Ranger ARC files and
    caches it to ``data/processed/`` (gitignored); later calls with the same
    ``(n_cells, seed)`` load straight from cache unless ``force_reload=True``.
    """
    cache_file = _cache_path(n_cells, seed)
    if cache_file.exists() and not force_reload:
        logger.info("Cache hit: loading %s", cache_file)
        return read_h5mu(cache_file)

    logger.info("Cache miss (%s): building from raw files", cache_file)
    mdata = _load_raw()
    mdata = subsample_mudata(mdata, n_cells=n_cells, seed=seed)

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    write_h5mu(cache_file, mdata)
    logger.info("Wrote cache to %s", cache_file)
    return mdata
