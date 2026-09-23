"""Test for the analysis-menu tool wrapper (uses the agent's small fixed-core cache)."""

from multiome_agent.menu.tf_motif_correlation import tf_motif_correlation


def test_tf_motif_correlation_known_tf(agent_fixed_core_mdata):
    result = tf_motif_correlation(agent_fixed_core_mdata, "SPI1")
    assert "error" not in result
    assert result["gene"] == "SPI1"
    assert -1.0 <= result["spearman_rho"] <= 1.0
    assert 0.0 <= result["pvalue"] <= 1.0


def test_tf_motif_correlation_unknown_gene_returns_error_not_exception(agent_fixed_core_mdata):
    result = tf_motif_correlation(agent_fixed_core_mdata, "NOT_A_REAL_GENE_XYZ")
    assert "error" in result
