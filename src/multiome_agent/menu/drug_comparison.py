"""Analysis-menu tools: treated-vs-control differential expression/
accessibility and batch-confound check for a drug/condition axis.

Thin agent-callable wrappers around the already-validated
`core.drug_comparison` functions -- per CLAUDE.md, the agent chooses *which*
validated analysis to run and with what parameters, it does not write new
analysis code.
"""

from __future__ import annotations

from multiome_agent.core.drug_comparison import (
    atac_differential_accessibility,
    batch_confound_check,
    rna_differential_expression,
)
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def drug_differential_expression(mdata, column_key: str, treated: str, control: str, cell_line: str | None = None) -> dict:
    """RNA differential expression between `treated` and `control` groups
    of `column_key` (e.g. "rna.drug_treatment", from
    `check_for_condition_groups`), optionally restricted to one real cell
    line's own cells. Returns an error dict (not a raised exception) for a
    missing column or an empty arm.
    """
    result = rna_differential_expression(mdata, column_key, treated, control, cell_line=cell_line)
    if "error" in result:
        logger.info("drug_differential_expression: %s", result["error"])
    else:
        logger.info(
            "drug_differential_expression: %s vs %s (cell_line=%s) -> %d significant genes",
            treated, control, cell_line, result["n_significant"],
        )
    return result


def drug_differential_accessibility(mdata, column_key: str, treated: str, control: str, cell_line: str | None = None) -> dict:
    """ATAC-side analog of `drug_differential_expression`: per-PEAK
    differential accessibility between the same two condition groups.
    """
    result = atac_differential_accessibility(mdata, column_key, treated, control, cell_line=cell_line)
    if "error" in result:
        logger.info("drug_differential_accessibility: %s", result["error"])
    else:
        logger.info(
            "drug_differential_accessibility: %s vs %s (cell_line=%s) -> %d significant peaks",
            treated, control, cell_line, result["n_significant"],
        )
    return result


def check_batch_confound(mdata, column_key: str, batch_column: str, cell_line: str | None = None) -> dict:
    """Whether `column_key`'s treated/control grouping is statistically
    confounded with `batch_column` (chi-squared test of independence on
    the real contingency table). Use this before trusting ANY treated-vs-
    control differential result -- if batch and treatment aren't
    independent, the result could be explained by batch alone.
    """
    result = batch_confound_check(mdata, column_key, batch_column, cell_line=cell_line)
    if "error" in result:
        logger.info("check_batch_confound: %s", result["error"])
    else:
        logger.info(
            "check_batch_confound: %s vs %s (cell_line=%s) -> likely_confound=%s (p=%s)",
            column_key, batch_column, cell_line, result["likely_confound"], result.get("pvalue"),
        )
    return result
