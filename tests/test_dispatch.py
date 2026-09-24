"""Tests for `data/dispatch.py`'s `load_fixed_core_from_local_config` --
the generic own-data entry point that auto-detects which
`config/local_paths.yaml` slot is populated, so a new user never has to
know or pass this project's internal 'tenx-cell-ranger' /
'shareseq-multi-cell-lines' dispatch-key strings for their own data.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

from multiome_agent.data.dispatch import load_fixed_core_from_local_config


def test_auto_detects_shareseq_config_when_both_rna_and_atac_paths_set():
    with patch("multiome_agent.data.dispatch.SHARESEQ_RNA_H5AD", "/fake/rna.h5ad"), \
         patch("multiome_agent.data.dispatch.SHARESEQ_ATAC_H5AD", "/fake/atac.h5ad"), \
         patch("multiome_agent.data.dispatch.TENX_MATRIX_H5", None), \
         patch("multiome_agent.data.dispatch.load_fixed_core_via_agent_decision") as mock_load:
        mock_load.return_value = ("mdata", "decision")
        result = load_fixed_core_from_local_config(model="claude-haiku-4-5")
        mock_load.assert_called_once_with("shareseq-multi-cell-lines", model="claude-haiku-4-5")
        assert result == ("mdata", "decision")


def test_auto_detects_tenx_config_when_only_combined_file_set():
    with patch("multiome_agent.data.dispatch.SHARESEQ_RNA_H5AD", None), \
         patch("multiome_agent.data.dispatch.SHARESEQ_ATAC_H5AD", None), \
         patch("multiome_agent.data.dispatch.TENX_MATRIX_H5", "/fake/combined.h5"), \
         patch("multiome_agent.data.dispatch.load_fixed_core_via_agent_decision") as mock_load:
        mock_load.return_value = ("mdata", "decision")
        result = load_fixed_core_from_local_config(model="claude-haiku-4-5")
        mock_load.assert_called_once_with("tenx-cell-ranger", model="claude-haiku-4-5")
        assert result == ("mdata", "decision")


def test_shareseq_config_takes_priority_when_both_are_set():
    with patch("multiome_agent.data.dispatch.SHARESEQ_RNA_H5AD", "/fake/rna.h5ad"), \
         patch("multiome_agent.data.dispatch.SHARESEQ_ATAC_H5AD", "/fake/atac.h5ad"), \
         patch("multiome_agent.data.dispatch.TENX_MATRIX_H5", "/fake/combined.h5"), \
         patch("multiome_agent.data.dispatch.load_fixed_core_via_agent_decision") as mock_load:
        mock_load.return_value = ("mdata", "decision")
        load_fixed_core_from_local_config()
        mock_load.assert_called_once_with("shareseq-multi-cell-lines", model=None)


def test_raises_clear_error_when_nothing_configured():
    with patch("multiome_agent.data.dispatch.SHARESEQ_RNA_H5AD", None), \
         patch("multiome_agent.data.dispatch.SHARESEQ_ATAC_H5AD", None), \
         patch("multiome_agent.data.dispatch.TENX_MATRIX_H5", None):
        with pytest.raises(ValueError, match="No local data configured"):
            load_fixed_core_from_local_config()
