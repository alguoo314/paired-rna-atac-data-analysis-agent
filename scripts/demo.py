"""One-command demo (`make demo`): loads the already-cached fixed-core
result and runs one real agent investigation end-to-end -- seconds and a
few cents, not the ~20-minute fixed-core pipeline itself (see README).

Ask your own question about the cached PBMC data:
    make demo QUESTION="does CD3E mark T cells in this dataset?"
    python scripts/demo.py "does CD3E mark T cells in this dataset?"
With no question given, falls back to a default SPI1 TF-motif-tracking question.
"""

from __future__ import annotations

import sys

from multiome_agent.agent.loop import run_agent

DEFAULT_QUESTION = (
    "In this PBMC dataset, does SPI1 expression correlate with its own "
    "motif's accessibility, and is that consistent with SPI1's known biology?"
)


def main() -> None:
    question = " ".join(sys.argv[1:]).strip() or DEFAULT_QUESTION
    print(f"Question: {question}\n")
    result = run_agent(question)
    print(result.answer)
    print(
        f"\n[turns={result.turn_count} cost=${result.estimated_cost_usd:.5f} "
        f"trail={result.reasoning_trail_path}]"
    )
    print(
        "\nThis is one live example. For the full picture -- fault-injection "
        "results, a literature-grounded known-biology checklist, and adversarially-"
        "judged novel findings -- see reports/examples/tenx-cell-ranger-report.md "
        "and shareseq-multi-cell-lines-report.md."
    )


if __name__ == "__main__":
    main()
