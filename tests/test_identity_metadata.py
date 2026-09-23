"""Tests for identity-metadata detection (core/identity_metadata.py). Tenx-cell-ranger
PBMC data has no such column (synthetic test only); shareseq-multi-cell-lines data really
does (cell_line_name) -- tested against the real cached fixed-core result.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData

from multiome_agent.core.identity_metadata import detect_identity_columns


def _make_mdata(obs: pd.DataFrame) -> MuData:
    n = len(obs)
    rna = AnnData(X=np.zeros((n, 3)), obs=obs)
    atac = AnnData(X=np.zeros((n, 3)), obs=obs.copy())
    return MuData({"rna": rna, "atac": atac})


def test_no_identity_column_reports_none_found():
    obs = pd.DataFrame({"leiden_rna": ["0", "1", "0", "1"], "sample": ["s1"] * 4})
    result = detect_identity_columns(_make_mdata(obs))
    assert result["none_found"] is True
    assert result["identity_columns_found"] == {}


def test_cell_line_name_column_detected_with_real_values():
    obs = pd.DataFrame({"cell_line_name": ["LineA"] * 3 + ["LineB"] * 2})
    result = detect_identity_columns(_make_mdata(obs))
    assert result["none_found"] is False
    assert result["identity_columns_found"]["rna.cell_line_name"] == {"LineA": 3, "LineB": 2}


def test_lineage_and_genotype_columns_also_detected():
    obs = pd.DataFrame({
        "atac_lineage": ["neural"] * 2 + ["epithelial"] * 3,
        "donor_id": ["D1", "D1", "D2", "D2", "D2"],
    })
    result = detect_identity_columns(_make_mdata(obs))
    assert "rna.atac_lineage" in result["identity_columns_found"]
    assert "rna.donor_id" in result["identity_columns_found"]


def test_depmap_column_detected_by_name():
    obs = pd.DataFrame({"Depmap": ["ACH-000001", "ACH-000001", "ACH-000958"]})
    result = detect_identity_columns(_make_mdata(obs))
    assert "rna.Depmap" in result["identity_columns_found"]
    assert result["identity_columns_found"]["rna.Depmap"] == {"ACH-000001": 2, "ACH-000958": 1}


def test_ach_id_values_detected_regardless_of_misleading_column_name():
    # Real bug caught against this project's actual shareseq-multi-cell-lines data: a column
    # literally named "Depmap" held unrelated batch labels ("b1"/"b2"/"b3"),
    # while the REAL DepMap ACH IDs lived in differently-named columns
    # ("rna_label"/"atac_label") that no name keyword would catch.
    obs = pd.DataFrame({
        "Depmap": ["b1", "b1", "b2"],
        "rna_label": ["ACH-000022", "ACH-000026", "ACH-000022"],
        "totally_unrelated_column": ["foo", "bar", "baz"],
    })
    result = detect_identity_columns(_make_mdata(obs))
    found = result["identity_columns_found"]
    assert "rna.rna_label" in found
    assert found["rna.rna_label"] == {"ACH-000022": 2, "ACH-000026": 1}
    assert "rna.Depmap" in found  # still surfaced by name, even though values aren't ACH IDs
    assert "rna.totally_unrelated_column" not in found


def test_real_shareseq_data_has_cell_line_name_with_real_counts():
    import pytest

    from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD

    if not (SHARESEQ_RNA_H5AD and SHARESEQ_ATAC_H5AD and SHARESEQ_RNA_HVG_H5AD):
        pytest.skip("config/local_paths.yaml not set up locally; shareseq-multi-cell-lines data unavailable")

    from multiome_agent.agent.shareseq_fixed_core_cache import get_shareseq_fixed_core

    mdata = get_shareseq_fixed_core()
    result = detect_identity_columns(mdata)
    assert result["none_found"] is False
    key = "rna.cell_line_name"
    assert key in result["identity_columns_found"]
    counts = result["identity_columns_found"][key]
    assert len(counts) == 8  # 8 real cell lines in this panel
    assert sum(counts.values()) == mdata.n_obs


def test_real_pbmc_data_has_no_identity_column(agent_fixed_core_mdata):
    result = detect_identity_columns(agent_fixed_core_mdata)
    assert result["none_found"] is True
