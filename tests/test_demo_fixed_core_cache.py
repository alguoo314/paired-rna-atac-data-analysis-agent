"""Tests for the small, committed demo cache (agent/demo_fixed_core_cache.py).

`get_demo_fixed_core` reads the already-committed `demo_cache/agent_demo_core.h5mu`
directly -- these tests don't rebuild it (that pays a real ~17-20 min cost,
see `scripts/build_demo_cache.py`), they verify the committed artifact is
still usable by every tool that reads it, and still free of the path leak
that motivated trimming it in the first place.
"""

from __future__ import annotations

from multiome_agent.agent.demo_fixed_core_cache import DEMO_N_CELLS, get_demo_fixed_core
from multiome_agent.core.clustering import ATAC_CLUSTER_KEY, RNA_CLUSTER_KEY
from multiome_agent.core.cross_modal_validation import cross_modal_marker_validation
from multiome_agent.menu.tf_motif_correlation import tf_motif_correlation


def test_demo_cache_loads_with_expected_shape():
    mdata = get_demo_fixed_core()
    assert mdata.n_obs == DEMO_N_CELLS
    assert RNA_CLUSTER_KEY in mdata.mod["rna"].obs
    assert ATAC_CLUSTER_KEY in mdata.mod["atac"].obs


def test_demo_cache_has_every_field_the_agent_tools_read():
    mdata = get_demo_fixed_core()
    rna, atac = mdata.mod["rna"], mdata.mod["atac"]
    assert "rank_genes_groups" in rna.uns
    assert "gene_activity" in atac.obsm
    assert "chromvar_deviations" in atac.obsm
    assert "chromvar_motif_names" in atac.uns
    assert "rna_atac_cluster_agreement" in mdata.uns


def test_demo_cache_tf_motif_correlation_still_works():
    mdata = get_demo_fixed_core()
    result = tf_motif_correlation(mdata, "SPI1")
    assert "error" not in result
    assert result["spearman_rho"] > 0.3  # a real, strong signal at this scale, not noise


def test_demo_cache_cross_modal_marker_check_still_works():
    mdata = get_demo_fixed_core()
    result = cross_modal_marker_validation(mdata, "CD3E")
    assert result["is_rna_marker"] is True


def test_demo_cache_has_no_leaked_local_file_paths():
    # Real bug this guards against: atac.uns["files"] once baked in this
    # machine's absolute fragments-file path -- caught before committing,
    # not hypothetical. Scan every uns dict for any path-like string.
    mdata = get_demo_fixed_core()
    for mod in ("rna", "atac"):
        for value in mdata.mod[mod].uns.values():
            text = str(value)
            assert "/ewsc/" not in text and "/home/" not in text and "/local_home/" not in text
