"""Drug/condition-axis detection: neither current dataset (tenx-cell-ranger, the
shareseq-multi-cell-lines panel) has a drug/treatment condition -- both are
control-only. This is a forward-looking branch for a future dataset that
does, per CLAUDE.md: the agent must determine control-only-vs-multi-
condition from the data itself, not be told, and if conditions exist, QC
on per-arm cell counts becomes relevant (a near-empty arm undermines any
condition-level claim).

Detects by COLUMN NAME keyword, not by inspecting/guessing at values --
deliberately narrow so this can't accidentally surface an identity-bearing
column (e.g. the shareseq-multi-cell-lines data's `cell_line_name`) through a generic
metadata dump; a drug/treatment column's VALUES (e.g. "DMSO" vs. "Drug_X")
aren't identity-sensitive the way cell-line names are, so reporting them
directly here is fine once such a column is actually found.
"""

from __future__ import annotations

from mudata import MuData

CONDITION_COLUMN_KEYWORDS = (
    "drug", "treatment", "condition", "dose", "compound", "perturbation", "guide", "sgrna", "stim",
)


def detect_condition_groups(mdata: MuData) -> dict:
    """Returns {"condition_columns_found": {"<modality>.<column>": {value: count}}, "is_control_only": bool}."""
    found: dict[str, dict[str, int]] = {}
    for mod_name in ("rna", "atac"):
        if mod_name not in mdata.mod:
            continue
        obs = mdata.mod[mod_name].obs
        for col in obs.columns:
            if any(kw in col.lower() for kw in CONDITION_COLUMN_KEYWORDS):
                counts = obs[col].value_counts()
                found[f"{mod_name}.{col}"] = {str(k): int(v) for k, v in counts.items()}
    return {"condition_columns_found": found, "is_control_only": len(found) == 0}


def condition_group_qc(mdata: MuData, column_key: str, min_cells: int = 50) -> dict:
    """Per-arm cell counts for one detected condition column (`column_key`
    as returned by `detect_condition_groups`, e.g. "rna.drug_treatment"),
    flagging any arm below `min_cells` -- a real concern the analysis
    should surface rather than silently averaging over an underpowered arm.
    """
    mod_name, col = column_key.split(".", 1)
    counts = mdata.mod[mod_name].obs[col].value_counts()
    underpowered = {str(k): int(v) for k, v in counts.items() if v < min_cells}
    return {
        "group_sizes": {str(k): int(v) for k, v in counts.items()},
        "underpowered_groups": underpowered,
        "any_underpowered": len(underpowered) > 0,
    }
