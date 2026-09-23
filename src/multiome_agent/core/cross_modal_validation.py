"""Cross-modal cell-type-call validation: the textbook ArchR/Signac use of
gene activity scores (see CLAUDE.md's design principles) -- given a marker
gene that's significant in one RNA cluster, does its ATAC gene-activity
independently confirm the SAME identity in that cluster's real cross-modal
partner (found via cell overlap, not by clusters coincidentally sharing an
integer label)? Two independently-derived signals agreeing is much stronger
identity evidence than either alone, and is the mechanism both the
known-biology checklist (section 2) and cell-identity discovery (section 3)
lean on.
"""

from __future__ import annotations

import anndata as ad
import pandas as pd
import scanpy as sc
from anndata import AnnData
from mudata import MuData

from multiome_agent.core.clustering import ATAC_CLUSTER_KEY, RNA_CLUSTER_KEY, gene_is_significant_marker_of_cluster


def best_matching_atac_cluster(mdata: MuData, rna_cluster: str) -> str:
    """The ATAC cluster containing the most cells from `rna_cluster` -- the
    natural cross-modal partner for validation, from cell overlap rather
    than clusters coincidentally sharing the same integer label (RNA and
    ATAC Leiden runs are independent, so "cluster 0" in each has no
    inherent relationship).
    """
    rna_labels = mdata.mod["rna"].obs[RNA_CLUSTER_KEY]
    atac_labels = mdata.mod["atac"].obs[ATAC_CLUSTER_KEY].reindex(rna_labels.index)
    contingency = pd.crosstab(rna_labels, atac_labels)
    return contingency.loc[rna_cluster].idxmax()


def gene_activity_is_significant_marker_of_cluster(
    atac: AnnData, gene: str, cluster: str, min_lfc: float = 0.5, max_padj: float = 0.05
) -> bool:
    """ATAC-side analog of `clustering.gene_is_significant_marker_of_cluster`,
    testing gene-activity accessibility instead of RNA expression. Lower
    default logFC threshold (0.5 vs. RNA's 1.0): gene activity is a noisier,
    indirect proxy (peak-window aggregation, not direct transcript counts),
    so real, biologically elevated signal typically shows a smaller effect
    size even when genuinely significant -- observed directly on real
    cross-modal marker checks, not assumed.
    """
    if gene not in atac.uns["gene_activity_genes"]:
        return False
    tmp = ad.AnnData(X=atac.obsm["gene_activity"], obs=atac.obs[[ATAC_CLUSTER_KEY]].copy())
    tmp.var_names = pd.Index(atac.uns["gene_activity_genes"])
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    sc.tl.rank_genes_groups(tmp, groupby=ATAC_CLUSTER_KEY, groups=[cluster], method="wilcoxon")
    df = sc.get.rank_genes_groups_df(tmp, group=cluster)
    row = df[df["names"] == gene]
    return bool(len(row) and row["pvals_adj"].iloc[0] < max_padj and row["logfoldchanges"].iloc[0] > min_lfc)


def cross_modal_marker_validation(mdata: MuData, gene: str) -> dict:
    """Is `gene` a significant RNA marker of some cluster, AND does its ATAC
    gene-activity independently confirm elevated accessibility in that
    cluster's real cross-modal partner? Returns enough detail to cite in a
    report even when the answer is "no" (which RNA cluster it marked,
    which ATAC cluster was checked) -- a validation attempt that finds no
    confirmation is still a real, reportable result, not a non-event.
    """
    rna = mdata.mod["rna"]
    marker_clusters = [
        c for c in rna.obs[RNA_CLUSTER_KEY].cat.categories
        if gene_is_significant_marker_of_cluster(rna, gene, c)
    ]
    if not marker_clusters:
        return {
            "gene": gene, "is_rna_marker": False, "rna_cluster": None,
            "matched_atac_cluster": None, "gene_activity_confirms": None,
        }

    rna_cluster = marker_clusters[0]
    atac_cluster = best_matching_atac_cluster(mdata, rna_cluster)
    confirms = gene_activity_is_significant_marker_of_cluster(mdata.mod["atac"], gene, atac_cluster)
    return {
        "gene": gene, "is_rna_marker": True, "rna_cluster": rna_cluster,
        "matched_atac_cluster": atac_cluster, "gene_activity_confirms": confirms,
    }
