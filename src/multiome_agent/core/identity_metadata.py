"""Identity-encoding metadata detection: does this dataset's own `.obs`
already carry a column that encodes cell-line/lineage/donor identity (e.g.
the shareseq-multi-cell-lines dataset's real, genotype-confirmed `cell_line_name`)?

This is legitimate discovery, not being told the answer: the agent finds it
by looking at what the data itself provides (column names, then real
values), the same way it would notice a `condition`/`drug` column via
`condition_detection.py`. It's a different, faster path to the SAME
identity-discovery goal as marker-based inference (`top_cluster_markers`,
`cross_modal_marker_check`) -- not a shortcut around it.

Two detection strategies, not just one:
1. COLUMN NAME keyword match (cheap, no values read until matched) --
   catches an explicitly-named column like `cell_line_name`.
2. COLUMN VALUE pattern match against a small sample -- catches an
   indirect identity encoding under an unhelpful or misleading name.
   Checked directly against this project's real shareseq-multi-cell-lines data: a column
   literally named "Depmap" turned out to hold unrelated batch labels
   ("b1"/"b2"/"b3"), while the REAL Broad Institute DepMap Model IDs
   ("ACH-XXXXXX") live in differently-named columns (`rna_label`,
   `atac_label`) that no name-keyword would have caught. Name-only
   detection would have missed exactly the column this was built for.

Deliberately no special-cased resolution here: this module only surfaces
matching columns and their raw values; recognizing an ACH-XXXXXX value as a
resolvable DepMap ID is left to the agent's own reasoning plus
`tools.depmap.resolve_depmap_id` (a real lookup against EBI's Cellosaurus
database), not hardcoded in this project's code.
"""

from __future__ import annotations

import re

from mudata import MuData

IDENTITY_COLUMN_KEYWORDS = (
    "cell_line", "cellline", "line_name", "lineage", "genotype", "donor_id", "patient_id", "cell_line_name",
    "depmap", "model_id",
)

# Broad Institute DepMap Model ID format -- the one concrete value-pattern
# this project has verified appears in real data under a non-obvious column
# name. Not an exhaustive "any ID-like string" heuristic (that would flag
# far too much); scoped to this one specific, real, checkable pattern.
_ACH_ID_RE = re.compile(r"^ACH-\d+$")
_VALUE_SAMPLE_SIZE = 20


def _column_values_look_like_depmap_ids(series) -> bool:
    try:
        sample = series.dropna().astype(str)
    except (TypeError, ValueError):
        return False
    if sample.empty:
        return False
    sample = sample.iloc[:_VALUE_SAMPLE_SIZE]
    return bool(len(sample)) and all(_ACH_ID_RE.match(v) for v in sample)


def detect_identity_columns(mdata: MuData) -> dict:
    """Returns {"identity_columns_found": {"<modality>.<column>": {value: count}}, "none_found": bool}.
    Real values are included (not anonymized) -- once the agent has found
    this itself from the data, per the project's cell-line-identity policy,
    using it is legitimate.
    """
    found: dict[str, dict[str, int]] = {}
    for mod_name in ("rna", "rna_hvg", "atac"):
        if mod_name not in mdata.mod:
            continue
        obs = mdata.mod[mod_name].obs
        for col in obs.columns:
            key = f"{mod_name}.{col}"
            name_match = any(kw in col.lower() for kw in IDENTITY_COLUMN_KEYWORDS)
            if not name_match and not _column_values_look_like_depmap_ids(obs[col]):
                continue
            counts = obs[col].value_counts()
            found[key] = {str(k): int(v) for k, v in counts.items()}
    return {"identity_columns_found": found, "none_found": len(found) == 0}
