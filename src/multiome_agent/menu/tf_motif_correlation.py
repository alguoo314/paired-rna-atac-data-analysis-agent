"""Analysis-menu tool: TF expression vs. its own motif's chromVAR deviation.

Thin agent-callable wrapper around the already-validated
`core.motif_deviations.tf_expression_motif_correlation` -- per CLAUDE.md,
the agent chooses *which* validated analysis to run and with what
parameters, it does not write new analysis code.
"""

from __future__ import annotations

from multiome_agent.core.motif_deviations import tf_expression_motif_correlation
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def tf_motif_correlation(mdata, gene: str) -> dict:
    """Spearman correlation between `gene`'s RNA expression and its own
    motif's chromVAR deviation score across cells. Returns an error dict
    (not a raised exception) if `gene` isn't in the data or has no motif
    match, so the agent loop can report it back to the model as a tool
    result instead of crashing.
    """
    gene = gene.strip().upper()
    rna = mdata["rna"]
    if gene not in rna.var_names:
        logger.info("tf_motif_correlation: gene %r not found in RNA data", gene)
        return {"error": f"Gene {gene!r} not found in this dataset's RNA var_names."}

    try:
        result = tf_expression_motif_correlation(mdata, gene, motif_name_contains=gene)
    except ValueError as e:
        logger.info("tf_motif_correlation: %s", e)
        return {"error": str(e)}

    logger.info(
        "tf_motif_correlation: %s vs %s -> rho=%.3f, p=%.2e",
        result["gene"], result["motif"], result["spearman_rho"], result["pvalue"],
    )
    return result
