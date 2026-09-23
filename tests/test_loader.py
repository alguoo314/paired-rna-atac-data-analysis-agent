"""Tests for the PBMC multiome data loader.

These exercise the real downloaded 10x dataset (not synthetic fixtures) since
the loader's job is specifically to parse that dataset's on-disk format
correctly. They skip gracefully if the raw files aren't present locally.
"""

import shutil

import pytest

from multiome_agent.data.loader import (
    CACHE_DIR,
    MATRIX_H5,
    PER_BARCODE_METRICS_CSV,
    PEAK_ANNOTATION_TSV,
    _cache_path,
    load_pbmc_multiome,
)

RAW_DATA_PRESENT = MATRIX_H5.exists() and PER_BARCODE_METRICS_CSV.exists() and PEAK_ANNOTATION_TSV.exists()

pytestmark = pytest.mark.skipif(
    not RAW_DATA_PRESENT,
    reason="raw PBMC multiome files not present in data/raw/pbmc10k_multiome/",
)

N_CELLS = 200


@pytest.fixture(scope="module")
def mdata():
    cache_file = _cache_path(N_CELLS, seed=0)
    if cache_file.exists():
        cache_file.unlink()
    return load_pbmc_multiome(n_cells=N_CELLS, seed=0, force_reload=True)


def test_rna_atac_split_has_right_feature_types(mdata):
    rna, atac = mdata["rna"], mdata["atac"]
    assert rna.n_vars > 0 and atac.n_vars > 0
    assert (rna.var["feature_types"] == "Gene Expression").all()
    assert (atac.var["feature_types"] == "Peaks").all()


def test_peak_coords_parsed(mdata):
    atac_var = mdata["atac"].var
    assert {"chrom", "start", "end"}.issubset(atac_var.columns)
    assert (atac_var["end"] > atac_var["start"]).all()
    # Most peaks sit on main chromosomes ("chr1", ...); GRCh38 also has
    # unplaced/alt scaffolds (e.g. "GL000194.1", "KI270713.1") without a
    # "chr" prefix, so we only assert the majority case, not all of them.
    assert atac_var["chrom"].str.startswith("chr").mean() > 0.9


def test_barcodes_match_between_modalities(mdata):
    assert list(mdata["rna"].obs_names) == list(mdata["atac"].obs_names)


def test_barcode_metrics_attached(mdata):
    for col in ["gex_umis_count", "atac_fragments", "atac_TSS_fragments", "atac_peak_region_fragments"]:
        assert col in mdata["rna"].obs.columns
        assert mdata["rna"].obs[col].notna().all()


def test_peak_annotation_attached(mdata):
    var = mdata["atac"].var
    assert "primary_gene" in var.columns
    assert var["primary_gene"].notna().any()


def test_subsample_returns_requested_cell_count(mdata):
    assert mdata.n_obs == N_CELLS


def test_cache_round_trip(tmp_path, monkeypatch):
    import multiome_agent.data.loader as loader_mod

    monkeypatch.setattr(loader_mod, "CACHE_DIR", tmp_path)

    n = 50
    cache_file = tmp_path / f"pbmc10k_multiome_n{n}_seed1.h5mu"
    assert not cache_file.exists()

    first = load_pbmc_multiome(n_cells=n, seed=1, force_reload=True)
    assert cache_file.exists()
    assert first.n_obs == n

    # Second call must hit the cache, not re-parse the raw h5: prove it by
    # making a raw reload raise, then confirming the cached call still succeeds.
    def _boom():
        raise AssertionError("cache miss: _load_raw() was called when it shouldn't have been")

    monkeypatch.setattr(loader_mod, "_load_raw", _boom)
    second = load_pbmc_multiome(n_cells=n, seed=1, force_reload=False)
    assert second.n_obs == n
    assert list(second.obs_names) == list(first.obs_names)
