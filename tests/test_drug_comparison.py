"""Tests for the drug/condition treated-vs-control comparison scaffold
(core/drug_comparison.py). Neither current dataset actually has a
drug/treatment condition axis -- this exercises the branch with
synthetically injected columns, same "inert scaffold + test" convention as
test_condition_detection.py.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData

from multiome_agent.core.drug_comparison import (
    atac_differential_accessibility,
    batch_confound_check,
    rna_differential_expression,
)


def _make_simple_drug_mdata():
    """20 control + 20 treated cells, no cell-line axis. UP_GENE is
    planted as a real, large treated-vs-control difference; FLAT_GENE has
    none (negative control). UP_PEAK mirrors UP_GENE on the ATAC side.
    """
    n_control, n_treated = 20, 20
    n = n_control + n_treated
    condition = ["control"] * n_control + ["treated"] * n_treated
    rng = np.random.RandomState(0)

    up_gene = np.concatenate([rng.poisson(2, n_control), rng.poisson(40, n_treated)]).astype(float)
    flat_gene = rng.poisson(10, n).astype(float)
    rna_X = np.column_stack([up_gene, flat_gene])
    rna_obs = pd.DataFrame({"drug_treatment": condition}, index=[f"cell{i}" for i in range(n)])
    rna = AnnData(X=rna_X, obs=rna_obs, var=pd.DataFrame(index=["UP_GENE", "FLAT_GENE"]))

    up_peak = np.concatenate([rng.poisson(1, n_control), rng.poisson(30, n_treated)]).astype(float)
    flat_peak = rng.poisson(8, n).astype(float)
    atac_X = np.column_stack([up_peak, flat_peak])
    atac = AnnData(X=atac_X, obs=rna_obs.copy(), var=pd.DataFrame(index=["UP_PEAK", "FLAT_PEAK"]))
    atac.layers["counts"] = atac_X.copy()

    return MuData({"rna": rna, "atac": atac})


def test_rna_differential_expression_detects_real_planted_difference():
    mdata = _make_simple_drug_mdata()
    result = rna_differential_expression(mdata, "rna.drug_treatment", "treated", "control")
    assert "error" not in result
    assert result["n_treated"] == 20 and result["n_control"] == 20
    up_genes = {r["gene"] for r in result["up_in_treated"]}
    assert "UP_GENE" in up_genes
    assert "FLAT_GENE" not in up_genes and "FLAT_GENE" not in {r["gene"] for r in result["down_in_treated"]}


def test_atac_differential_accessibility_detects_real_planted_difference():
    mdata = _make_simple_drug_mdata()
    result = atac_differential_accessibility(mdata, "rna.drug_treatment", "treated", "control")
    assert "error" not in result
    up_peaks = {r["peak"] for r in result["up_in_treated"]}
    assert "UP_PEAK" in up_peaks
    assert "FLAT_PEAK" not in up_peaks


def test_rna_differential_expression_missing_column_returns_error():
    mdata = _make_simple_drug_mdata()
    result = rna_differential_expression(mdata, "rna.not_a_real_column", "treated", "control")
    assert "error" in result


def test_rna_differential_expression_empty_group_returns_error():
    mdata = _make_simple_drug_mdata()
    result = rna_differential_expression(mdata, "rna.drug_treatment", "treated", "NOT_A_REAL_GROUP")
    assert "error" in result


def _make_cell_line_confound_drug_mdata():
    """The real failure mode pooling risks: MARKER_A is purely a cell-line
    marker (high in LINE_A, low in LINE_B, with NO real condition effect
    within either line) -- but the treated arm is composition-skewed
    toward LINE_A and the control arm toward LINE_B, so a POOLED
    treated-vs-control test sees MARKER_A as significantly "up in
    treated" purely from cell-line composition, not a real drug effect.
    Scoped to one line alone, the confound must disappear.
    """
    # LINE_A: 15 treated + 5 control, MARKER_A always high (~20, small jitter).
    # LINE_B: 5 treated + 15 control, MARKER_A always low (~2, small jitter).
    rng = np.random.RandomState(0)
    n_a_treated, n_a_control, n_b_treated, n_b_control = 15, 5, 5, 15
    cell_line = (["LINE_A"] * (n_a_treated + n_a_control)) + (["LINE_B"] * (n_b_treated + n_b_control))
    condition = (
        ["treated"] * n_a_treated + ["control"] * n_a_control
        + ["treated"] * n_b_treated + ["control"] * n_b_control
    )
    n = len(cell_line)
    marker_a = np.concatenate([
        rng.poisson(20, n_a_treated + n_a_control),  # LINE_A: always high, regardless of condition
        rng.poisson(2, n_b_treated + n_b_control),   # LINE_B: always low, regardless of condition
    ]).astype(float)

    rna_obs = pd.DataFrame(
        {"drug_treatment": condition, "cell_line_name": cell_line}, index=[f"cell{i}" for i in range(n)],
    )
    rna = AnnData(X=marker_a.reshape(-1, 1), obs=rna_obs, var=pd.DataFrame(index=["MARKER_A"]))
    atac = AnnData(X=np.zeros((n, 1)), obs=rna_obs.copy(), var=pd.DataFrame(index=["PEAK1"]))
    atac.layers["counts"] = atac.X.copy()
    return MuData({"rna": rna, "atac": atac})


def test_rna_differential_expression_pooled_finds_a_spurious_cell_line_confound():
    mdata = _make_cell_line_confound_drug_mdata()
    result = rna_differential_expression(mdata, "rna.drug_treatment", "treated", "control")
    assert "error" not in result
    assert "MARKER_A" in {r["gene"] for r in result["up_in_treated"]}


def test_rna_differential_expression_cell_line_scoping_excludes_the_confound():
    mdata = _make_cell_line_confound_drug_mdata()
    for line in ("LINE_A", "LINE_B"):
        result = rna_differential_expression(mdata, "rna.drug_treatment", "treated", "control", cell_line=line)
        assert "error" not in result
        all_flagged = {r["gene"] for r in result["up_in_treated"]} | {r["gene"] for r in result["down_in_treated"]}
        assert "MARKER_A" not in all_flagged


def test_rna_differential_expression_cell_line_no_such_column_returns_error():
    mdata = _make_simple_drug_mdata()  # has no cell_line_name column at all
    result = rna_differential_expression(mdata, "rna.drug_treatment", "treated", "control", cell_line="LINE_A")
    assert "error" in result


def test_batch_confound_check_detects_real_confound():
    # Batch perfectly predicts treatment arm -- the real failure mode a
    # Judger reviewing a drug-comparison finding needs to catch.
    n = 40
    condition = ["control"] * 20 + ["treated"] * 20
    batch = ["batch1"] * 20 + ["batch2"] * 20
    rna_obs = pd.DataFrame({"drug_treatment": condition, "batch": batch}, index=[f"cell{i}" for i in range(n)])
    rna = AnnData(X=np.zeros((n, 1)), obs=rna_obs)
    atac = AnnData(X=np.zeros((n, 1)), obs=rna_obs.copy())
    mdata = MuData({"rna": rna, "atac": atac})

    result = batch_confound_check(mdata, "rna.drug_treatment", "batch")
    assert "error" not in result
    assert result["likely_confound"] is True
    assert result["pvalue"] < 0.05


def test_batch_confound_check_no_confound_when_batch_is_balanced():
    n = 40
    condition = (["control", "treated"] * 20)
    batch = (["batch1"] * 20 + ["batch2"] * 20)
    rna_obs = pd.DataFrame({"drug_treatment": condition, "batch": batch}, index=[f"cell{i}" for i in range(n)])
    rna = AnnData(X=np.zeros((n, 1)), obs=rna_obs)
    atac = AnnData(X=np.zeros((n, 1)), obs=rna_obs.copy())
    mdata = MuData({"rna": rna, "atac": atac})

    result = batch_confound_check(mdata, "rna.drug_treatment", "batch")
    assert "error" not in result
    assert result["likely_confound"] is False


def test_batch_confound_check_missing_batch_column_returns_error():
    n = 10
    rna_obs = pd.DataFrame({"drug_treatment": ["control"] * 5 + ["treated"] * 5}, index=[f"cell{i}" for i in range(n)])
    rna = AnnData(X=np.zeros((n, 1)), obs=rna_obs)
    atac = AnnData(X=np.zeros((n, 1)), obs=rna_obs.copy())
    mdata = MuData({"rna": rna, "atac": atac})

    result = batch_confound_check(mdata, "rna.drug_treatment", "not_a_real_batch_column")
    assert "error" in result
