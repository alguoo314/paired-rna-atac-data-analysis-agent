from pathlib import Path

from multiome_agent import config
from multiome_agent.data.loader import _resolve_matrix_h5


def test_agent_model_default():
    assert config.AGENT_MODEL == "claude-haiku-4-5"


def test_max_turns_default():
    assert config.MAX_TURNS == 15
    assert isinstance(config.MAX_TURNS, int)


def test_real_env_api_key_loads():
    # The project .env is expected to have a real key set during development.
    assert isinstance(config.ANTHROPIC_API_KEY, str)
    assert len(config.ANTHROPIC_API_KEY) > 10


def test_missing_local_paths_file_returns_none(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "_LOCAL_PATHS_FILE", tmp_path / "does_not_exist.yaml")
    assert config._load_local_paths() == {}


def test_shareseq_paths_load_when_local_config_present_else_are_none():
    # `test_missing_local_paths_file_returns_none` above covers the loading
    # function's behavior when config/local_paths.yaml is absent, in
    # isolation. This dev environment has that file configured (the shareseq-multi-cell-lines
    # multi-cell-line data phase), so the real module-level constants are
    # non-None strings here -- asserting they're always None was only ever
    # true by coincidence of an earlier dev environment's state, not a
    # property of the code. Assert the shape (None, or a non-empty string),
    # not a hardcoded presence/absence that depends on local machine state.
    for path in (config.SHARESEQ_RNA_H5AD, config.SHARESEQ_ATAC_H5AD, config.SHARESEQ_RNA_HVG_H5AD):
        assert path is None or (isinstance(path, str) and len(path) > 0)


def test_tenx_matrix_h5_override_loads_from_local_paths_yaml(tmp_path, monkeypatch):
    yaml_path = tmp_path / "local_paths.yaml"
    yaml_path.write_text("tenx_matrix_h5: /some/custom/path.h5\n")
    monkeypatch.setattr(config, "_LOCAL_PATHS_FILE", yaml_path)
    assert config._load_local_paths()["tenx_matrix_h5"] == "/some/custom/path.h5"


def test_tenx_matrix_h5_is_none_shape_when_unset():
    # Same "shape, not hardcoded presence" reasoning as the shareseq-multi-cell-lines paths
    # test above -- this dev environment may or may not have the key set.
    assert config.TENX_MATRIX_H5 is None or (isinstance(config.TENX_MATRIX_H5, str) and len(config.TENX_MATRIX_H5) > 0)


def test_resolve_matrix_h5_falls_back_to_raw_dir_when_no_override():
    result = _resolve_matrix_h5(None, Path("/some/raw/dir"))
    assert result == Path("/some/raw/dir/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")


def test_resolve_matrix_h5_uses_override_when_set():
    # Real fix for a real limitation (previously an unconditionally
    # hardcoded path, flagged in the README): a user's own combined 10x
    # .h5 file, anywhere, via config/local_paths.yaml's tenx_matrix_h5.
    result = _resolve_matrix_h5("/my/own/data/custom_name.h5", Path("/some/raw/dir"))
    assert result == Path("/my/own/data/custom_name.h5")
