"""Lets the agent inspect raw data file structure and decide which loader a
given data source needs, instead of a human/script pre-selecting "combined"
vs. "two-file" loading logic (see `agent.loader_selector`). Reads only cheap
metadata via `h5py` directly -- shapes, column names, small categorical
value sets -- never the matrix data itself, so nothing here violates
CLAUDE.md's "summaries, not raw data, go to the LLM" principle.

Two real on-disk formats exist in this project (verified directly against
both, not assumed):
- 10x Cell Ranger ARC's own HDF5 layout (`/matrix/...`, used by the tenx-cell-ranger
  PBMC file): one file holds BOTH modalities, distinguished by
  `matrix/features/feature_type` containing both "Gene Expression" and
  "Peaks".
- Standard AnnData `.h5ad` layout (`/X`, `/obs`, `/var`, used by the shareseq-multi-cell-lines
  data's files): each file is single-modality, no `feature_types` column at
  all in this project's shareseq files.
"""

from __future__ import annotations

import h5py


def _inspect_one(path: str) -> dict:
    with h5py.File(path, "r") as f:
        top_level = list(f.keys())

        if "matrix" in top_level:
            shape = f["matrix/shape"][:]  # 10x convention: [n_vars, n_obs]
            raw_types = f["matrix/features/feature_type"][:]
            feature_types = sorted({t.decode() if isinstance(t, bytes) else t for t in raw_types})
            return {
                "file_format": "10x_cellranger_h5",
                "n_vars": int(shape[0]),
                "n_obs": int(shape[1]),
                "feature_types_present": feature_types,
            }

        if {"X", "obs", "var"}.issubset(top_level):
            x_shape = f["X"].attrs.get("shape")
            var_columns = list(f["var"].keys())
            obs_columns = list(f["obs"].keys())
            return {
                "file_format": "anndata_h5ad",
                "n_obs": int(x_shape[0]) if x_shape is not None else None,
                "n_vars": int(x_shape[1]) if x_shape is not None else None,
                "var_columns": var_columns,
                "obs_columns": obs_columns,
                "has_feature_types_column": "feature_types" in var_columns,
            }

        return {"file_format": "unknown", "top_level_keys": top_level}


def inspect_data_paths(paths: list[str]) -> dict:
    """Structural summary of each given file -- format, shape, and (for
    combined-format files) which modalities its features span. Raises if a
    path doesn't exist or isn't a readable HDF5 file, rather than silently
    returning a partial report.
    """
    return {path: _inspect_one(path) for path in paths}
