"""Fixed-core clustering: RNA Leiden + markers, ATAC LSI + Leiden, cluster agreement.

CLAUDE.md: "Discordance is expected, not a bug" -- RNA and ATAC correlate
weakly per gene, so cluster agreement (ARI) is measured and reported, not
optimized for. A low ARI is a normal finding, not a bug to chase.
"""

from __future__ import annotations

import anndata as ad
import muon as mu
import pandas as pd
import scanpy as sc
from anndata import AnnData
from mudata import MuData
from sklearn.metrics import adjusted_rand_score

from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

RNA_CLUSTER_KEY = "leiden_rna"
ATAC_CLUSTER_KEY = "leiden_atac"


def cluster_rna(rna: AnnData, n_top_genes: int = 2000, resolution: float = 1.0) -> AnnData:
    """Standard scanpy pipeline: normalize -> HVG -> PCA -> neighbors -> Leiden -> markers.

    Runs on a copy of counts in `.layers["counts"]` so `.X` ends up
    log-normalized (scanpy's expected state for `rank_genes_groups`), while
    the raw counts stay available for anything downstream that needs them.
    """
    sc.pp.filter_genes(rna, min_cells=3)
    rna.layers["counts"] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=n_top_genes)
    sc.pp.pca(rna, mask_var="highly_variable")
    sc.pp.neighbors(rna)
    sc.tl.leiden(rna, resolution=resolution, key_added=RNA_CLUSTER_KEY, flavor="igraph")
    sc.tl.rank_genes_groups(rna, groupby=RNA_CLUSTER_KEY, method="wilcoxon")

    logger.info(
        "RNA clustering: %d cells -> %d Leiden clusters",
        rna.n_obs, rna.obs[RNA_CLUSTER_KEY].nunique(),
    )
    return rna


def cluster_atac(atac: AnnData, resolution: float = 1.0) -> AnnData:
    """TF-IDF -> LSI (muon) -> neighbors -> Leiden, the standard scATAC pipeline."""
    # Peaks with zero counts in this subsample would divide-by-zero in IDF.
    sc.pp.filter_genes(atac, min_cells=3)
    atac.layers["counts"] = atac.X.copy()
    mu.atac.pp.tfidf(atac)
    mu.atac.tl.lsi(atac)
    # First LSI component usually tracks sequencing depth, not biology (standard
    # scATAC practice -- e.g. Signac/ArchR both drop it before downstream steps).
    # `n_pcs` in sc.pp.neighbors takes the *first* N columns of `use_rep`, so
    # dropping component 0 needs an explicit slice, not a smaller `n_pcs`.
    atac.obsm["X_lsi_no_depth"] = atac.obsm["X_lsi"][:, 1:]
    sc.pp.neighbors(atac, use_rep="X_lsi_no_depth")
    sc.tl.leiden(atac, resolution=resolution, key_added=ATAC_CLUSTER_KEY, flavor="igraph")

    logger.info(
        "ATAC clustering: %d cells -> %d Leiden clusters",
        atac.n_obs, atac.obs[ATAC_CLUSTER_KEY].nunique(),
    )
    return atac


def gene_is_significant_marker_of_cluster(
    rna: AnnData, gene: str, cluster: str, min_lfc: float = 1.0, max_padj: float = 0.05
) -> bool:
    """Whether `gene` is a significant positive marker of this one `cluster`.

    Checks `pvals_adj`/`logfoldchanges` from `rank_genes_groups` rather than
    literal top-N rank: with real data, a canonical marker (e.g. CD14) can be
    genuinely, significantly enriched in the right cluster while still
    sitting outside an arbitrary top-25 cutoff by raw score -- especially
    when Leiden resolution splits one cell type across a few clusters. The
    padj/logFC thresholds are the standard "is this a marker" criterion
    (matches Seurat/Scanpy convention), so they're more robust to clustering
    granularity than a rank cutoff.
    """
    df = sc.get.rank_genes_groups_df(rna, group=cluster)
    row = df[df["names"] == gene]
    return bool(len(row) and row["pvals_adj"].iloc[0] < max_padj and row["logfoldchanges"].iloc[0] > min_lfc)


def marker_is_recovered(
    rna: AnnData, gene: str, min_lfc: float = 1.0, max_padj: float = 0.05
) -> bool:
    """Whether `gene` is a significant positive marker of at least one cluster."""
    return any(
        gene_is_significant_marker_of_cluster(rna, gene, cluster, min_lfc, max_padj)
        for cluster in rna.obs[RNA_CLUSTER_KEY].cat.categories
    )


def top_rna_cluster_markers(rna: AnnData, cluster: str, n: int = 15) -> list[dict]:
    """Top `n` positive markers of `cluster` by `rank_genes_groups` score --
    open-ended discovery (which genes distinguish this cluster) rather than
    `gene_is_significant_marker_of_cluster`'s yes/no check against one
    pre-specified gene. This is the real evidence an agent needs to infer
    what a cluster IS (cell type / cell line / condition) from the data
    itself instead of being told directly -- see agent/prompts.py.
    """
    df = sc.get.rank_genes_groups_df(rna, group=cluster).sort_values("scores", ascending=False).head(n)
    return df[["names", "logfoldchanges", "pvals_adj"]].rename(columns={"names": "gene"}).to_dict("records")


def top_gene_activity_cluster_markers(atac: AnnData, cluster: str, n: int = 15) -> list[dict]:
    """Same idea as `top_rna_cluster_markers`, computed on ATAC gene-activity
    scores (see core.gene_activity) instead of RNA expression -- the
    ATAC-side, independently-derived analog. Builds a small throwaway
    AnnData over `gene_activity` (already sparse, so this stays memory-light
    even at ATAC's larger gene counts) and reuses `rank_genes_groups`, the
    same validated method the RNA side uses, rather than a hand-rolled
    differential test.
    """
    tmp = ad.AnnData(X=atac.obsm["gene_activity"], obs=atac.obs[[ATAC_CLUSTER_KEY]].copy())
    tmp.var_names = pd.Index(atac.uns["gene_activity_genes"])
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    sc.tl.rank_genes_groups(tmp, groupby=ATAC_CLUSTER_KEY, groups=[cluster], method="wilcoxon")
    df = sc.get.rank_genes_groups_df(tmp, group=cluster).sort_values("scores", ascending=False).head(n)
    return df[["names", "logfoldchanges", "pvals_adj"]].rename(columns={"names": "gene"}).to_dict("records")


def rna_atac_cluster_agreement(mdata: MuData) -> dict:
    """Adjusted Rand index + contingency table between RNA and ATAC Leiden labels.

    Returns a dict with keys `"ari"` and `"contingency"` rather than just
    logging, so a later report-writer step has a concrete table to cite.
    """
    rna_labels = mdata["rna"].obs[RNA_CLUSTER_KEY]
    atac_labels = mdata["atac"].obs[ATAC_CLUSTER_KEY].reindex(rna_labels.index)

    ari = adjusted_rand_score(rna_labels, atac_labels)
    contingency = pd.crosstab(rna_labels, atac_labels, rownames=["rna"], colnames=["atac"])

    logger.info("RNA-vs-ATAC cluster agreement: ARI=%.3f (see contingency table for detail)", ari)
    return {"ari": ari, "contingency": contingency}
