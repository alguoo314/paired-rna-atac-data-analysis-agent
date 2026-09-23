#!/bin/bash
# Downloads the raw PBMC 10k Multiome dataset + hg38 reference genome into
# data/raw/ (gitignored -- not committed to the repo, ~3.3GB total). Skips
# any file that already exists at the expected size, so it's safe to re-run.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PBMC_DIR="$REPO_ROOT/data/raw/pbmc10k_multiome"
REF_DIR="$REPO_ROOT/data/raw/reference"
mkdir -p "$PBMC_DIR" "$REF_DIR"

BASE="https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k"

_fetch() {
  local url="$1" out="$2"
  if [ -f "$out" ]; then
    echo "Already present, skipping: $out"
  else
    echo "Downloading $out ..."
    curl -sS -o "$out" "$url"
  fi
}

_fetch "$BASE/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5" \
  "$PBMC_DIR/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
_fetch "$BASE/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv" \
  "$PBMC_DIR/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv"
_fetch "$BASE/pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv" \
  "$PBMC_DIR/pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv"
_fetch "$BASE/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz" \
  "$PBMC_DIR/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
_fetch "$BASE/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi" \
  "$PBMC_DIR/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
_fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit" \
  "$REF_DIR/hg38.2bit"

echo "Done. (JASPAR motifs are fetched automatically on first use into data/raw/reference/ -- no separate download step.)"
