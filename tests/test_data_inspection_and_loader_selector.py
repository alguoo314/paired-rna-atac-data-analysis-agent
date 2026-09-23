"""Tests for the agent-mediated loader-selection path (see
`multiome_agent.tools.data_inspection` and `multiome_agent.agent.loader_selector`).
The structural-inspection tests are pure/offline; the decision test hits the
real Anthropic API on Haiku (a single small call) and is skipped without a
key, matching this repo's existing convention (see test_agent_loop.py).
"""

from __future__ import annotations

import os

import pytest

from multiome_agent.config import SHARESEQ_ATAC_H5AD, SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD
from multiome_agent.data.loader import MATRIX_H5
from multiome_agent.tools.data_inspection import inspect_data_paths

PUBLIC_RAW_AVAILABLE = MATRIX_H5.exists()
SHARESEQ_DATA_CONFIGURED = bool(SHARESEQ_RNA_H5AD and SHARESEQ_ATAC_H5AD and SHARESEQ_RNA_HVG_H5AD)


@pytest.mark.skipif(not PUBLIC_RAW_AVAILABLE, reason="raw PBMC 10x h5 not downloaded locally")
def test_inspect_public_combined_file_reveals_both_feature_types():
    report = inspect_data_paths([str(MATRIX_H5)])
    info = report[str(MATRIX_H5)]
    assert info["file_format"] == "10x_cellranger_h5"
    assert set(info["feature_types_present"]) == {"Gene Expression", "Peaks"}
    assert info["n_obs"] > 0 and info["n_vars"] > 0


@pytest.mark.skipif(not SHARESEQ_DATA_CONFIGURED, reason="config/local_paths.yaml not set up locally")
def test_inspect_shareseq_two_files_have_no_combined_feature_types_column():
    report = inspect_data_paths([SHARESEQ_RNA_H5AD, SHARESEQ_ATAC_H5AD])
    for path in (SHARESEQ_RNA_H5AD, SHARESEQ_ATAC_H5AD):
        info = report[path]
        assert info["file_format"] == "anndata_h5ad"
        assert info["has_feature_types_column"] is False
        assert info["n_obs"] == 5814  # both files share the same cell count (same-cell multiome)


@pytest.mark.skipif(not os.environ.get("ANTHROPIC_API_KEY"), reason="no ANTHROPIC_API_KEY available")
@pytest.mark.skipif(not PUBLIC_RAW_AVAILABLE, reason="raw PBMC 10x h5 not downloaded locally")
def test_agent_correctly_decides_combined_strategy_for_public_file():
    from multiome_agent.agent.loader_selector import decide_loading_strategy

    decision = decide_loading_strategy([str(MATRIX_H5)], model="claude-haiku-4-5")
    assert decision.strategy == "combined_single_file"
    assert decision.cost_usd > 0


@pytest.mark.skipif(not os.environ.get("ANTHROPIC_API_KEY"), reason="no ANTHROPIC_API_KEY available")
@pytest.mark.skipif(not SHARESEQ_DATA_CONFIGURED, reason="config/local_paths.yaml not set up locally")
def test_agent_correctly_decides_two_file_strategy_for_shareseq_files():
    from multiome_agent.agent.loader_selector import decide_loading_strategy

    decision = decide_loading_strategy(
        [SHARESEQ_RNA_H5AD, SHARESEQ_RNA_HVG_H5AD, SHARESEQ_ATAC_H5AD], model="claude-haiku-4-5"
    )
    assert decision.strategy == "separate_per_modality_files"
    assert decision.cost_usd > 0
