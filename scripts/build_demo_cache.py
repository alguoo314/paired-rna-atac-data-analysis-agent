"""Regenerates `demo_cache/agent_demo_core.h5mu` (see
`agent/demo_fixed_core_cache.py` for what it keeps/drops and why). Not run
by end users -- `make demo` loads the already-committed result directly.
Only needed if this cache is ever intentionally rebuilt (e.g. after a
pipeline change). Real cost: pays the ~17-20 min fixed fragments-import
cost once (see `core/pipeline.py`).
"""

from multiome_agent.agent.demo_fixed_core_cache import build_demo_cache


def main() -> None:
    build_demo_cache()
    print("Done.")


if __name__ == "__main__":
    main()
