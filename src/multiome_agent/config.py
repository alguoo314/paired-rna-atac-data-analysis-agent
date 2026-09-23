"""Central config: env vars (.env) + optional local paths (config/local_paths.yaml)."""

import os
from pathlib import Path

import yaml
from dotenv import load_dotenv

REPO_ROOT = Path(__file__).resolve().parents[2]

load_dotenv(REPO_ROOT / ".env")

AGENT_MODEL = os.environ.get("AGENT_MODEL", "claude-haiku-4-5")
MAX_TURNS = int(os.environ.get("AGENT_MAX_TURNS", "15"))
LOG_LEVEL = os.environ.get("LOG_LEVEL", "INFO")


def get_anthropic_api_key() -> str:
    """Read the key lazily so importing config never fails just because it's unset."""
    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        raise RuntimeError(
            "ANTHROPIC_API_KEY is not set. Add it to the project .env file."
        )
    return key


# Kept as a module attribute too, for callers that only need to check presence/display,
# not exercise the lazy error above.
ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY")

_LOCAL_PATHS_FILE = REPO_ROOT / "config" / "local_paths.yaml"


def _load_local_paths() -> dict:
    if not _LOCAL_PATHS_FILE.exists():
        return {}
    with open(_LOCAL_PATHS_FILE) as f:
        return yaml.safe_load(f) or {}


_local_paths = _load_local_paths()

# Absolute paths to the owner's shareseq-multi-cell-lines multiome data (see CLAUDE.md).
# None until config/local_paths.yaml is created locally; that file is gitignored
# because it names filesystem locations, not because it holds secrets.
SHARESEQ_RNA_H5AD = _local_paths.get("shareseq_rna_h5ad")
SHARESEQ_ATAC_H5AD = _local_paths.get("shareseq_atac_h5ad")
SHARESEQ_RNA_HVG_H5AD = _local_paths.get("shareseq_rna_hvg_h5ad")

# Optional override for the tenx-cell-ranger combined-single-file 10x .h5 path
# -- None (the default) means `data/loader.py` falls back to the standard
# `download_pbmc_data.sh` location. Lets a user point this at their own
# combined 10x Cell Ranger ARC .h5 file without editing code, the same way
# the shareseq per-modality paths above are configured.
TENX_MATRIX_H5 = _local_paths.get("tenx_matrix_h5")
