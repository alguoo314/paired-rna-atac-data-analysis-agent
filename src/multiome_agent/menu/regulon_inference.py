"""Analysis-menu tool: regulon inference (TF -> predicted target-gene-set).

Thin agent-callable wrapper around the already-validated
`core.regulon_inference.infer_regulon_targets` -- per CLAUDE.md, the agent
chooses *which* validated analysis to run and with what TF, it does not
write new analysis code.
"""

from __future__ import annotations

from multiome_agent.core.regulon_inference import infer_regulon_targets
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def regulon_inference(mdata, tf_gene: str, target_gene: str | None = None, cell_line: str | None = None) -> dict:
    """Find candidate target genes for `tf_gene` via real motif + chromatin
    evidence (motif present in a target's own promoter, or in a distal peak
    that itself significantly links to that target's expression), then
    report, PER TARGET GENE, the Spearman correlation between the TF's own
    RNA and that target's own RNA -- never an aggregate score, never the
    TF's own expression standing in for a target. Returns an error dict
    (not a raised exception) if the TF isn't in the RNA data, has no
    matching motif, or `motif_match` data isn't available in this dataset
    (compute_motif_deviations must have run and persisted it).

    `targets` is capped (see `infer_regulon_targets`'s docstring for why --
    an abundant motif can have thousands of real candidates, too many to
    return in full). If you already have a SPECIFIC target gene in mind
    (e.g. checking a literature-reported TF-target pair), pass
    `target_gene` -- its real entry is guaranteed to appear in `targets`
    even if it wouldn't otherwise make the cap.

    `cell_line`, if given, restricts the correlation to only that cell
    line's own cells (see `core.regulon_inference.infer_regulon_targets`'s
    docstring for why pooling across a multi-cell-line dataset is the
    wrong default -- it can report a spurious between-line confound as a
    real TF-target relationship).

    `tf_detection_rate` (top-level) and each target's own
    `target_detection_rate` report the fraction of cells with nonzero
    expression, scoped the same way as the correlation -- ALWAYS check
    these before treating a significant rho/q as real: a gene detected in
    only a handful of cells can't support a meaningful correlation
    regardless of how significant it looks, and this is a DIFFERENT, more
    decisive check than cluster-marker status (see
    `infer_regulon_targets`'s docstring).
    """
    result = infer_regulon_targets(mdata, tf_gene, target_gene=target_gene, cell_line=cell_line)
    if "error" in result:
        logger.info("regulon_inference: %s", result["error"])
    else:
        logger.info(
            "regulon_inference: %s -> %d candidate targets, any_significant=%s",
            result["tf"], result["n_candidate_targets"], result["any_significant_target"],
        )
    return result
