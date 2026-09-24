"""Scores the fixed-core pipeline's output against the known-biology
checklist (`eval/checklists/pbmc_known_biology.yaml`, written as ground
truth before any of this scoring code existed -- see PROGRESS.md step 7).

Every recovered/not-recovered verdict carries the actual computed numbers
behind it (rho, p-value, which cluster, etc.), not just a bare boolean --
CLAUDE.md wants every claim grounded in a specific computed table.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc
import yaml
from mudata import MuData
from scipy.stats import mannwhitneyu

from multiome_agent.config import REPO_ROOT
from multiome_agent.core.clustering import RNA_CLUSTER_KEY, gene_is_significant_marker_of_cluster
from multiome_agent.core.motif_deviations import tf_expression_motif_correlation
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

CHECKLIST_PATH = REPO_ROOT / "src" / "multiome_agent" / "eval" / "checklists" / "pbmc_known_biology.yaml"


def load_checklist(path: Path = CHECKLIST_PATH) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def annotate_clusters(rna, checklist: dict) -> dict[str, str | None]:
    """Label each RNA Leiden cluster with a checklist cell_type, only when
    exactly one distinct checklist cell_type has a significant marker hit
    for it. Ambiguous (>1 cell_type hits) or unmatched clusters get `None`
    -- this dataset's clusters aren't all checklist cell types (no
    platelet/DC entries in the checklist), so leaving some unlabeled is
    expected, not an error.
    """
    labels: dict[str, str | None] = {}
    for cluster in rna.obs[RNA_CLUSTER_KEY].cat.categories:
        hit_cell_types = {
            marker["cell_type"]
            for marker in checklist["rna_markers"]
            if gene_is_significant_marker_of_cluster(rna, marker["gene"], cluster)
        }
        labels[cluster] = hit_cell_types.pop() if len(hit_cell_types) == 1 else None
    return labels


def score_rna_markers(rna, checklist: dict) -> list[dict]:
    results = []
    for marker in checklist["rna_markers"]:
        gene, cell_type = marker["gene"], marker["cell_type"]
        hit_clusters = [
            c for c in rna.obs[RNA_CLUSTER_KEY].cat.categories
            if gene_is_significant_marker_of_cluster(rna, gene, c)
        ]
        results.append({
            "gene": gene, "cell_type": cell_type,
            "recovered": len(hit_clusters) > 0,
            "significant_in_clusters": hit_clusters,
        })
    return results


def _motif_names_to_check(motif_item: dict) -> list[str]:
    names = [motif_item["motif_name"]]
    for alt in motif_item.get("motif_family_alternates", []):
        # "MA0466.4/CEBPB" -> "CEBPB"
        names.append(alt.split("/")[-1])
    return names


def score_motifs(mdata: MuData, checklist: dict, cluster_labels: dict[str, str | None]) -> list[dict]:
    atac = mdata.mod["atac"]
    deviations = atac.obsm["chromvar_deviations"]
    rna_cluster_of_cell = mdata.mod["rna"].obs[RNA_CLUSTER_KEY].reindex(deviations.index)

    results = []
    for motif_item in checklist["motifs"]:
        cell_type = motif_item["cell_type"]
        candidate_names = _motif_names_to_check(motif_item)
        matched_cols = [
            col for col in deviations.columns
            if any(name.lower() in col.lower() for name in candidate_names)
        ]

        entry = {
            "motif_id": motif_item["motif_id"], "motif_name": motif_item["motif_name"],
            "cell_type": cell_type, "matched_columns": matched_cols,
        }
        if not matched_cols:
            entry["status"] = "not_recovered"
            entry["reason"] = "no matching motif column in chromvar_deviations"
            results.append(entry)
            continue

        labeled_clusters = [c for c, lbl in cluster_labels.items() if lbl == cell_type]
        if not labeled_clusters:
            entry["status"] = "present_but_not_scoreable"
            entry["reason"] = f"no cluster was uniquely labeled '{cell_type}'"
            results.append(entry)
            continue

        # If multiple JASPAR columns matched (e.g. a family alternate also
        # substring-matches another motif's name), test each and report the
        # single best (lowest p-value) rather than silently picking one.
        best = None
        for col in matched_cols:
            in_cluster = deviations.loc[rna_cluster_of_cell.isin(labeled_clusters), col]
            rest = deviations.loc[~rna_cluster_of_cell.isin(labeled_clusters), col]
            stat, pval = mannwhitneyu(in_cluster, rest, alternative="greater")
            candidate = {
                "matched_column": col, "labeled_clusters": labeled_clusters,
                "pvalue": float(pval), "median_in_cluster": float(in_cluster.median()),
                "median_rest": float(rest.median()),
            }
            if best is None or candidate["pvalue"] < best["pvalue"]:
                best = candidate

        entry.update(best)
        entry["status"] = "recovered" if best["pvalue"] < 0.05 else "not_recovered"
        results.append(entry)
    return results


def score_tf_expression_tracks_motif(mdata: MuData, checklist: dict) -> list[dict]:
    results = []
    for item in checklist["tf_expression_tracks_motif"]:
        gene = item["gene"]
        try:
            corr = tf_expression_motif_correlation(mdata, gene=gene, motif_name_contains=item["motif_name"])
        except ValueError as e:
            results.append({"gene": gene, "recovered": False, "reason": str(e)})
            continue
        recovered = corr["spearman_rho"] > 0 and corr["pvalue"] < 0.05
        results.append({
            "gene": gene, "motif": corr["motif"], "spearman_rho": corr["spearman_rho"],
            "pvalue": corr["pvalue"], "recovered": recovered,
        })
    return results


def score_checklist(mdata: MuData, checklist_path: Path = CHECKLIST_PATH) -> dict:
    """Full scorecard: RNA markers, motifs, TF-expression-tracks-motif, plus
    recall summary counts. `mdata` must already have `run_fixed_core` applied.
    """
    checklist = load_checklist(checklist_path)
    rna = mdata.mod["rna"]

    cluster_labels = annotate_clusters(rna, checklist)
    rna_results = score_rna_markers(rna, checklist)
    motif_results = score_motifs(mdata, checklist, cluster_labels)
    tf_results = score_tf_expression_tracks_motif(mdata, checklist)

    def _recall_bool(results):
        n = len(results)
        hit = sum(1 for r in results if r["recovered"])
        return {"hit": hit, "total": n, "recall": hit / n if n else float("nan")}

    def _recall_status(results):
        n = len(results)
        hit = sum(1 for r in results if r["status"] == "recovered")
        return {"hit": hit, "total": n, "recall": hit / n if n else float("nan")}

    scorecard = {
        "cluster_labels": cluster_labels,
        "rna_markers": rna_results,
        "motifs": motif_results,
        "tf_expression_tracks_motif": tf_results,
        "summary": {
            "rna_markers": _recall_bool(rna_results),
            "motifs": _recall_status(motif_results),
            "tf_expression_tracks_motif": _recall_bool(tf_results),
        },
    }
    logger.info(
        "Checklist scorecard: RNA markers %d/%d, motifs %d/%d, TF-motif %d/%d",
        scorecard["summary"]["rna_markers"]["hit"], scorecard["summary"]["rna_markers"]["total"],
        scorecard["summary"]["motifs"]["hit"], scorecard["summary"]["motifs"]["total"],
        scorecard["summary"]["tf_expression_tracks_motif"]["hit"], scorecard["summary"]["tf_expression_tracks_motif"]["total"],
    )
    return scorecard
