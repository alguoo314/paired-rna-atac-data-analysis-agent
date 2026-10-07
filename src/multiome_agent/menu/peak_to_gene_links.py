"""Analysis-menu tool: peak-to-gene links (cis-co-accessibility).

Thin agent-callable wrapper around the already-validated
`core.peak_to_gene_links.peak_to_gene_links` -- per CLAUDE.md, the agent
chooses *which* validated analysis to run and with what gene, it does not
write new analysis code.
"""

from __future__ import annotations

from multiome_agent.core.peak_to_gene_links import peak_to_gene_links as _peak_to_gene_links
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def peak_to_gene_links(mdata, gene: str) -> dict:
    """Test whether any DISTAL peak's accessibility (i.e. NOT already
    captured by `gene`'s own gene-activity score) tracks `gene`'s own RNA
    expression across cells. Returns an error dict (not a raised exception)
    if `gene` isn't in the RNA data or has no protein-coding GENCODE
    coordinates, so the agent loop can report it back to the model as a
    tool result instead of crashing.
    """
    result = _peak_to_gene_links(mdata, gene)
    if "error" in result:
        logger.info("peak_to_gene_links: %s", result["error"])
    else:
        logger.info(
            "peak_to_gene_links: %s -> %d candidate distal peaks, any_significant=%s",
            result["gene"], result["n_candidate_distal_peaks"], result["any_significant_distal_link"],
        )
    return result
