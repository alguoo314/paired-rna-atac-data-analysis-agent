"""Stretch-phase results plot: cell-line-label-swap fault severity vs. the
pipeline's cell-line-recovery ARI, on the private multi-cell-line dataset.

Real numbers hardcoded from PROGRESS.md's "Stretch -- Private data:
cell-line/mixed-sample fault types" section (the fault-injection functions
that produced them need the private dataset locally to rerun, which most
readers of this repo won't have -- this script plots the already-verified,
committed numbers rather than requiring a rerun). Statistics only, no
per-cell or identity data.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_PATH = REPO_ROOT / "reports" / "examples" / "stretch_phase_results.png"

SEVERITY = [0.0, 0.10, 0.30, 0.50]
RECOVERY_ARI = [0.767, 0.623, 0.358, 0.157]


def main() -> None:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(SEVERITY, RECOVERY_ARI, marker="o", markersize=8, zorder=3)
    for x, y in zip(SEVERITY, RECOVERY_ARI):
        ax.annotate(f"{y:.3f}", (x, y), textcoords="offset points", xytext=(8, 6), fontsize=9)

    ax.set_xlabel("Fraction of cells with swapped cell-line label")
    ax.set_ylabel("Cell-line recovery ARI (RNA clustering vs. true label)")
    ax.set_ylim(0, 0.85)
    ax.set_title("Cell-line-label-swap fault: severity vs. recovery signal\n(private multi-cell-line dataset)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150)
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
