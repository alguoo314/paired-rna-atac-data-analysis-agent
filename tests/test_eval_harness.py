"""Cheap tests for the eval harness's classifiers and qc_summary composition
-- no real API calls, no fixed-core reruns. `build_scenarios` itself
(which includes a real ~17-min snapatac2 rerun) is exercised for real by
the comprehensive report (`agent/report_generator.py`), not re-run here.
"""

from __future__ import annotations

from multiome_agent.eval.eval_harness import (
    _classify_detected,
    _classify_diagnosis_matches,
    _qc_summary_for_scenario,
)


def test_classify_detected_true_on_problem_language():
    assert _classify_detected("There is a clear problem: the doublet rate looks elevated.")


def test_classify_detected_false_on_clean_language():
    assert not _classify_detected("This dataset looks reasonably healthy overall. No major red flags.")


def test_classify_detected_false_on_other_negated_phrasing():
    assert not _classify_detected("Nothing here is unusual, and there is no issue with quality.")


def test_classify_detected_true_when_a_real_problem_follows_clean_language():
    text = "No major issues overall, though the doublet rate is unusually elevated in a subset of cells."
    assert _classify_detected(text)


def test_classify_diagnosis_matches_clean_control_is_none():
    assert _classify_diagnosis_matches("anything", "clean_control") is None


def test_classify_diagnosis_matches_doublets_keyword():
    assert _classify_diagnosis_matches("The elevated doublet score suggests multiplets.", "injected_doublets")


def test_classify_diagnosis_matches_wrong_keyword_is_false():
    assert not _classify_diagnosis_matches("Fragment depth looks low.", "injected_doublets")


class _FakeMdata:
    class _Mod(dict):
        def __getitem__(self, key):
            return super().__getitem__(key)

    def __init__(self):
        import numpy as np
        import pandas as pd
        from anndata import AnnData

        rna_obs = pd.DataFrame({
            "n_genes_by_counts": [1800, 1900], "total_counts": [3700, 3800],
            "pct_counts_mt": [9.0, 10.0], "predicted_doublet": [False, True],
            "doublet_score": [0.02, 0.5], "leiden_rna": pd.Categorical(["0", "1"]),
        })
        atac_obs = pd.DataFrame({
            "frip": [0.7, 0.8], "tss_enrichment": [15.0, 17.0], "nucleosome_signal": [0.9, 0.95],
            "atac_fragments": [12000, 13000], "leiden_atac": pd.Categorical(["0", "1"]),
        })
        self.mod = {
            "rna": AnnData(X=np.zeros((2, 3)), obs=rna_obs),
            "atac": AnnData(X=np.zeros((2, 3)), obs=atac_obs),
        }
        self.uns = {"rna_atac_cluster_agreement": {"ari": 0.5}}
        self.n_obs = 2

    def __getitem__(self, key):
        return self.mod[key]


def _patch_spi1_rho(monkeypatch, rho=0.5, pval=1e-10):
    import multiome_agent.agent.qc_summary as qcs

    monkeypatch.setattr(
        qcs, "tf_expression_motif_correlation",
        lambda mdata, gene, motif_name_contains: {"gene": gene, "motif": "fake", "spearman_rho": rho, "pvalue": pval},
    )


def test_qc_summary_for_clean_scenario_has_no_additional_note(monkeypatch):
    _patch_spi1_rho(monkeypatch)
    mdata = _FakeMdata()
    text = _qc_summary_for_scenario(mdata, "clean_control", None)
    assert "Additional note" not in text
    assert "ARI" in text or "cluster agreement" in text.lower()


def test_qc_summary_for_doublets_scenario_is_one_blended_summary_with_no_subgroup_leak(monkeypatch):
    """Regression test for the same category of "cheating" bug as the
    baseline-leak fix: the doublets qc_summary must not name a "core"/
    "extra"/"real"/"synthetic" split anywhere -- just one blended aggregate
    over the full augmented population, exactly like an ordinary QC summary."""
    _patch_spi1_rho(monkeypatch)
    mdata = _FakeMdata()
    fault_score = {
        "blended_n_cells": 556,
        "blended_median_genes_per_cell": 1800.0, "blended_median_umis_per_cell": 3700.0,
        "blended_median_pct_mt": 9.5, "blended_predicted_doublet_rate": 0.15,
        "blended_median_doublet_score": 0.08,
    }
    text = _qc_summary_for_scenario(mdata, "injected_doublets", fault_score)
    assert "extra" not in text.lower()
    assert "core cells" not in text.lower() and "core cell" not in text.lower()
    assert "synthetic" not in text.lower()
    assert "556 cells" in text
    assert "15.0% predicted doublets" in text
    assert "0.080" in text


def test_qc_summary_for_shuffled_pairing_has_no_baseline_leak(monkeypatch):
    """Regression test for the "cheating" bug: the faulted summary must show
    ONLY the faulted ARI/rho, with no separate clean-baseline number
    anywhere in the text a model could use to infer "something changed"
    without its own judgment."""
    _patch_spi1_rho(monkeypatch, rho=0.9)  # clean value, must NOT appear below
    mdata = _FakeMdata()
    fault_score = {"ari_faulted": 0.123, "ari_clean": 0.9, "tf_motif_rho_faulted": 0.05, "tf_motif_rho_clean": 0.9}
    text = _qc_summary_for_scenario(mdata, "shuffled_rna_atac_pairing", fault_score)
    assert "baseline" not in text.lower()
    assert "0.123" in text
    assert "0.900" not in text  # the clean rho (0.9) must not leak in anywhere


def test_qc_summary_for_downsampling_has_no_baseline_leak(monkeypatch):
    _patch_spi1_rho(monkeypatch)
    mdata = _FakeMdata()
    fault_score = {
        "n_fragment_median_faulted": 2656.0, "n_fragment_median_clean": 13310.0,
        "tss_enrichment_median_faulted": 14.48, "tss_enrichment_median_clean": 16.7,
        "nucleosome_signal_median_faulted": 0.93, "nucleosome_signal_median_clean": 0.9,
    }
    text = _qc_summary_for_scenario(mdata, "atac_downsampling", fault_score)
    assert "baseline" not in text.lower()
    assert "2656" in text
    assert "13310" not in text and "13,310" not in text


def test_run_model_eval_passes_model_through_without_real_api_call(monkeypatch):
    """Confirms `run_model_eval`/`_run_scenario` thread `model` into
    `run_agent` and into the returned `ScenarioResult`, using a stub so this
    costs nothing."""
    import multiome_agent.eval.eval_harness as eh

    class _FakeResult:
        answer = "The data looks clean, no problems."
        turn_count = 1
        estimated_cost_usd = 0.001
        reasoning_trail_path = "/tmp/fake.md"
        tool_calls = []
        refused = False
        refusal_detail = None

    seen_models = []

    def _fake_run_agent(question, model=None, mdata=None, qc_summary=None, **kw):
        seen_models.append(model)
        return _FakeResult()

    monkeypatch.setattr(eh, "run_agent", _fake_run_agent)
    scenarios = [{"fault_type": "clean_control", "severity": 0.0, "description": "clean", "qc_summary": "x"}]
    results = eh.run_model_eval("claude-opus-5", scenarios, clean_mdata=object(), n_consistency_repeats=0)

    assert results["model"] == "claude-opus-5"
    assert seen_models == ["claude-opus-5"]
    assert results["scenario_results"][0]["refused"] is False
