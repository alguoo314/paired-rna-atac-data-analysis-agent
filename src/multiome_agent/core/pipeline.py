"""Fixed-core orchestrator: RNA QC -> ATAC QC -> clustering -> cluster agreement.

One entry point (`run_fixed_core`) for the always-run baseline pass that
CLAUDE.md's later agent loop and report writer build on. The individual
pieces (`rna_qc.run_rna_qc`, `atac_qc.run_atac_qc`, `clustering.cluster_rna`,
etc.) stay independently importable/composable for anything that wants only
one piece (e.g. tests, or a future step reusing just the clustering).
"""

from __future__ import annotations

from pathlib import Path

from mudata import MuData

from multiome_agent.core.atac_qc import run_atac_qc
from multiome_agent.core.clustering import cluster_atac, cluster_rna, rna_atac_cluster_agreement
from multiome_agent.core.gene_activity import compute_gene_activity
from multiome_agent.core.motif_deviations import compute_motif_deviations, tf_expression_motif_correlation
from multiome_agent.core.rna_qc import run_rna_qc
from multiome_agent.core.snap_import import import_snap_fragments
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)


def run_fixed_core(mdata: MuData, fragments_file: Path) -> MuData:
    """Run the fixed-core QC + clustering + gene-activity + motif pipeline in
    place; returns `mdata`.

    Cluster agreement (ARI + contingency table) is stashed in
    `mdata.uns["rna_atac_cluster_agreement"]` since it's a cross-modality
    result with no natural home in either modality's `.obs`/`.var`.

    The fragments file is imported via snapatac2 exactly once (~a few
    minutes, dominated by sorting the whole file) and reused for both TSS
    enrichment (ATAC QC) and gene activity, instead of importing it twice.
    """
    logger.info("Running fixed-core pipeline on %d cells", mdata.n_obs)

    run_rna_qc(mdata["rna"])

    snap_data = import_snap_fragments(mdata["atac"].obs_names, fragments_file)
    run_atac_qc(mdata["atac"], fragments_file, snap_data=snap_data)
    compute_gene_activity(mdata, snap_data)

    cluster_rna(mdata["rna"])
    cluster_atac(mdata["atac"])
    mdata.uns["rna_atac_cluster_agreement"] = rna_atac_cluster_agreement(mdata)

    compute_motif_deviations(mdata)
    # SPI1/PU.1 is a well-characterized PBMC myeloid TF (CLAUDE.md's known-biology
    # notes) with an unambiguous single-name JASPAR motif match -- one representative
    # sanity check, not an exhaustive TF audit.
    tf_expression_motif_correlation(mdata, gene="SPI1", motif_name_contains="Spi1")

    logger.info("Fixed-core pipeline complete")
    return mdata
