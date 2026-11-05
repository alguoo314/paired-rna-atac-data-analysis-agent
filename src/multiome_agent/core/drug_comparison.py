"""Treated-vs-control differential expression/accessibility, scoped per
cell line -- the "peak-level differential accessibility between two named
groups" item previously flagged in CLAUDE.md's roadmap as not yet
implemented, now real.

A forward-looking branch, same status as `core.condition_detection`: NEITHER
current dataset (tenx-cell-ranger, shareseq-multi-cell-lines) actually has a
drug/treatment condition axis -- both are control-only, confirmed by
`detect_condition_groups`, not assumed. This is exercised only via
synthetic tests today; a future dataset that does have a condition axis is
what this is built for.

Scoped per cell line for the SAME reason `peak_to_gene_links`/
`regulon_inference` are (see their docstrings, and PROGRESS_phase2.md):
pooling a treated-vs-control comparison across multiple distinct cell lines
risks confounding real treatment effect with between-line differences in
cell-line COMPOSITION across the two arms -- a gene that's simply a marker
of whichever line happens to be over-represented in the "treated" arm would
look like a real drug effect, pooled, with zero real within-line treatment
signal at all.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from mudata import MuData
from scipy.stats import chi2_contingency


def _condition_group_mask(mdata: MuData, column_key: str, treated_value: str, control_value: str, cell_line: str | None):
    """Shared setup for `rna_differential_expression`/`atac_differential_accessibility`:
    resolves the real group labels (reindexed onto RNA's own cell order,
    since the condition column can live on either modality), the combined
    treated-or-control mask, and the optional cell_line mask -- same
    `cell_line_name`-based scoping convention as `core.peak_to_gene_links`/
    `core.regulon_inference`'s own `cell_line` parameter.
    """
    rna = mdata.mod["rna"]
    mod_name, col = column_key.split(".", 1)
    if mod_name not in mdata.mod or col not in mdata.mod[mod_name].obs.columns:
        return None, None, f"Condition column {column_key!r} not found in this dataset."
    group_col = mdata.mod[mod_name].obs[col].reindex(rna.obs_names)

    mask = group_col.isin([treated_value, control_value]).to_numpy()
    if cell_line is not None:
        if "cell_line_name" not in rna.obs.columns:
            return None, None, "This dataset has no 'cell_line_name' column -- cell_line scoping isn't applicable here."
        mask = mask & (rna.obs["cell_line_name"] == cell_line).to_numpy()
    return group_col, mask, None


def _differential_result(adata, group_col, mask, treated_value, control_value, min_lfc, max_padj, name_col):
    n_treated = int((group_col[mask] == treated_value).sum())
    n_control = int((group_col[mask] == control_value).sum())
    if n_treated == 0 or n_control == 0:
        return {"error": f"No cells for one or both groups ({treated_value!r}, {control_value!r}) in this scope."}

    sub = adata[mask].copy()
    sub.obs["_drug_comparison_group"] = group_col[mask].astype(str).to_numpy()
    sc.tl.rank_genes_groups(
        sub, groupby="_drug_comparison_group", groups=[treated_value],
        reference=control_value, method="wilcoxon",
    )
    df = sc.get.rank_genes_groups_df(sub, group=treated_value)
    df["significant"] = (df["pvals_adj"] < max_padj) & (df["logfoldchanges"].abs() > min_lfc)
    sig = df[df["significant"]].copy()
    sig["_abs_lfc"] = sig["logfoldchanges"].abs()
    sig = sig.sort_values("_abs_lfc", ascending=False)

    def _records(sub_df):
        return (
            sub_df.head(15)[["names", "logfoldchanges", "pvals_adj"]]
            .rename(columns={"names": name_col})
            .to_dict("records")
        )

    return {
        "n_treated": n_treated, "n_control": n_control,
        "n_significant": int(sig.shape[0]),
        "up_in_treated": _records(sig[sig["logfoldchanges"] > 0]),
        "down_in_treated": _records(sig[sig["logfoldchanges"] < 0]),
    }


def rna_differential_expression(
    mdata: MuData, column_key: str, treated_value: str, control_value: str,
    cell_line: str | None = None, min_lfc: float = 0.5, max_padj: float = 0.05,
) -> dict:
    """RNA differential expression between `treated_value` and
    `control_value` groups of condition column `column_key` (e.g.
    "rna.drug_treatment", as returned by `detect_condition_groups`),
    optionally restricted to one real cell line's own cells first (see
    module docstring for why pooling across lines is the wrong default on
    a multi-cell-line dataset). Returns the top up/down genes by effect
    size, plus real counts for both arms, so a downstream literature-
    checking step gets specific, named genes to look up -- never an
    aggregate "N genes changed" summary.

    Returns an error dict (not a raised exception) for a missing column,
    missing cell_line scoping, or an empty arm -- same convention as
    `peak_to_gene_links`/`infer_regulon_targets`.
    """
    group_col, mask, error = _condition_group_mask(mdata, column_key, treated_value, control_value, cell_line)
    if error:
        return {"error": error}
    result = _differential_result(mdata.mod["rna"], group_col, mask, treated_value, control_value, min_lfc, max_padj, "gene")
    if "error" in result:
        return result
    result.update({"column_key": column_key, "treated": treated_value, "control": control_value, "cell_line": cell_line})
    return result


def atac_differential_accessibility(
    mdata: MuData, column_key: str, treated_value: str, control_value: str,
    cell_line: str | None = None, min_lfc: float = 0.5, max_padj: float = 0.05,
) -> dict:
    """ATAC-side analog of `rna_differential_expression`: per-PEAK
    differential accessibility (not gene activity -- real peak-level
    accessibility, the roadmap item this implements) between the same two
    condition groups, using `atac.layers["counts"]` (real per-cell peak
    counts) if present, falling back to `.X` otherwise -- same convention
    as `peak_to_gene_links`.
    """
    rna = mdata.mod["rna"]
    group_col, mask, error = _condition_group_mask(mdata, column_key, treated_value, control_value, cell_line)
    if error:
        return {"error": error}

    atac = mdata.mod["atac"]
    atac_for_test = atac.copy()
    if "counts" in atac_for_test.layers:
        atac_for_test.X = atac_for_test.layers["counts"]
    # The condition column's own index may be ATAC's or RNA's -- reindex the
    # mask (built against RNA's cell order in `_condition_group_mask`) onto
    # ATAC's own cell order, same barcode-alignment convention the rest of
    # this project's cross-modal code already relies on.
    atac_order_mask = pd.Series(mask, index=rna.obs_names).reindex(atac_for_test.obs_names).fillna(False).to_numpy()
    atac_group_col = group_col.reindex(atac_for_test.obs_names)

    result = _differential_result(
        atac_for_test, atac_group_col, atac_order_mask, treated_value, control_value, min_lfc, max_padj, "peak",
    )
    if "error" in result:
        return result
    result.update({"column_key": column_key, "treated": treated_value, "control": control_value, "cell_line": cell_line})
    return result


def batch_confound_check(mdata: MuData, column_key: str, batch_column: str, cell_line: str | None = None) -> dict:
    """Real, deterministic check for whether `column_key`'s treated/control
    grouping is confounded with `batch_column` (e.g. "rna.batch") -- i.e.
    whether batch and treatment assignment are statistically independent,
    via a chi-squared test of independence on the real contingency table.
    If they're NOT independent (small p-value), any treated-vs-control
    differential result on this scope could be explained by batch alone,
    not a real drug effect -- exactly what an adversarial Judger reviewing
    a drug-comparison finding needs to check independently (per CLAUDE.md
    principle 1: the LLM reasons about this REAL number, it never guesses
    or generates one).

    Returns an error dict (not a raised exception) for a missing column.
    """
    rna = mdata.mod["rna"]
    mod_name, col = column_key.split(".", 1)
    if mod_name not in mdata.mod or col not in mdata.mod[mod_name].obs.columns:
        return {"error": f"Condition column {column_key!r} not found in this dataset."}
    if batch_column not in rna.obs.columns:
        return {"error": f"Batch column {batch_column!r} not found in this dataset's RNA .obs."}

    condition = mdata.mod[mod_name].obs[col].reindex(rna.obs_names)
    batch = rna.obs[batch_column]
    mask = pd.Series(True, index=rna.obs_names)
    if cell_line is not None:
        if "cell_line_name" not in rna.obs.columns:
            return {"error": "This dataset has no 'cell_line_name' column -- cell_line scoping isn't applicable here."}
        mask = rna.obs["cell_line_name"] == cell_line

    table = pd.crosstab(condition[mask], batch[mask])
    if table.shape[0] < 2 or table.shape[1] < 2:
        return {
            "column_key": column_key, "batch_column": batch_column, "cell_line": cell_line,
            "contingency_table": table.to_dict(), "pvalue": None, "likely_confound": False,
            "note": "Fewer than 2 real groups on one axis -- confound test isn't meaningful here.",
        }

    chi2, pvalue, dof, _ = chi2_contingency(table)
    likely_confound = bool(pvalue < 0.05)
    return {
        "column_key": column_key, "batch_column": batch_column, "cell_line": cell_line,
        "contingency_table": table.to_dict(), "chi2": float(chi2), "pvalue": float(pvalue),
        "likely_confound": likely_confound,
    }
