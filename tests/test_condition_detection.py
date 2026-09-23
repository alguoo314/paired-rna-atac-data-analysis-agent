"""Tests for the drug/condition-axis detection scaffold (core/condition_detection.py).
Neither current dataset actually has a condition column -- this exercises
the branch with a synthetically injected one, per the "inert scaffold +
test" plan, plus confirms real data correctly reports control-only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData

from multiome_agent.core.condition_detection import condition_group_qc, detect_condition_groups


def _make_mdata(rna_obs: pd.DataFrame, atac_obs: pd.DataFrame | None = None) -> MuData:
    n = len(rna_obs)
    rna = AnnData(X=np.zeros((n, 3)), obs=rna_obs)
    atac = AnnData(X=np.zeros((n, 3)), obs=atac_obs if atac_obs is not None else rna_obs.copy())
    return MuData({"rna": rna, "atac": atac})


def test_no_condition_column_reports_control_only():
    obs = pd.DataFrame({"leiden_rna": ["0", "1", "0", "1"], "cell_line_name": ["A", "B", "A", "B"]})
    mdata = _make_mdata(obs)
    result = detect_condition_groups(mdata)
    assert result["is_control_only"] is True
    assert result["condition_columns_found"] == {}


def test_injected_drug_column_is_detected_with_real_group_sizes():
    obs = pd.DataFrame({
        "leiden_rna": ["0"] * 30 + ["1"] * 70,
        "drug_treatment": ["DMSO"] * 60 + ["Drug_X_10uM"] * 40,
    })
    mdata = _make_mdata(obs)
    result = detect_condition_groups(mdata)
    assert result["is_control_only"] is False
    assert "rna.drug_treatment" in result["condition_columns_found"]
    sizes = result["condition_columns_found"]["rna.drug_treatment"]
    assert sizes["DMSO"] == 60 and sizes["Drug_X_10uM"] == 40


def test_condition_group_qc_flags_underpowered_arm():
    obs = pd.DataFrame({"condition": ["control"] * 95 + ["treated"] * 5})
    mdata = _make_mdata(obs)
    result = condition_group_qc(mdata, "rna.condition", min_cells=50)
    assert result["group_sizes"] == {"control": 95, "treated": 5}
    assert result["any_underpowered"] is True
    assert result["underpowered_groups"] == {"treated": 5}


def test_condition_group_qc_no_flag_when_all_arms_large_enough():
    obs = pd.DataFrame({"condition": ["control"] * 60 + ["treated"] * 55})
    mdata = _make_mdata(obs)
    result = condition_group_qc(mdata, "rna.condition", min_cells=50)
    assert result["any_underpowered"] is False


def test_case_insensitive_and_multiple_keyword_matches():
    obs = pd.DataFrame({
        "leiden_rna": ["0"] * 10,
        "Compound_Dose": ["low"] * 5 + ["high"] * 5,
        "sgRNA_target": ["NT"] * 5 + ["GENE1"] * 5,
    })
    mdata = _make_mdata(obs)
    result = detect_condition_groups(mdata)
    assert "rna.Compound_Dose" in result["condition_columns_found"]
    assert "rna.sgRNA_target" in result["condition_columns_found"]


def test_real_pbmc_and_shareseq_data_are_control_only(agent_fixed_core_mdata):
    result = detect_condition_groups(agent_fixed_core_mdata)
    assert result["is_control_only"] is True
