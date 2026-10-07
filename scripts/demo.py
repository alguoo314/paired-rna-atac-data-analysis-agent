"""One-command demo (`make demo`): loads the small, committed demo cache
(`demo_cache/agent_demo_core.h5mu`, 150 cells, ~46MB -- see
`agent/demo_fixed_core_cache.py`) and runs one real agent investigation
end-to-end -- seconds and a few cents, genuinely on a fresh clone, not the
~17-20 minute full fixed-core pipeline (see README). This is illustrative
only, at a small subsample -- see `reports/examples/` for the real,
full-scale analyses this project's flagship reports are built from.

Ask your own question about the cached PBMC data:
    make demo QUESTION="does CD3E mark T cells in this dataset?"
    python scripts/demo.py "does CD3E mark T cells in this dataset?"
With no question given, falls back to a default GATA3 regulon-inference question --
deliberately showcasing `regulon_inference` (a TF tracking SPECIFIC other genes it's
predicted to regulate, not just its own motif) rather than the narrower `tf_motif_correlation`
check a plain TF-vs-own-motif question would exercise.
"""

from __future__ import annotations

import sys

from multiome_agent.agent.demo_fixed_core_cache import get_demo_fixed_core
from multiome_agent.agent.loop import run_agent
from multiome_agent.agent.qc_summary import fixed_core_summary, format_qc_summary

DEFAULT_QUESTION = (
    "Does GATA3 regulate any specific target gene in this dataset, not just its own motif?"
)


def main() -> None:
    question = " ".join(sys.argv[1:]).strip() or DEFAULT_QUESTION
    print(f"Question: {question}\n")
    mdata = get_demo_fixed_core()
    qc_summary = format_qc_summary(fixed_core_summary(mdata))
    result = run_agent(question, mdata=mdata, qc_summary=qc_summary)
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
