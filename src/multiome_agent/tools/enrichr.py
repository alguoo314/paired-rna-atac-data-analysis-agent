"""Enrichr gene-set enrichment, cached locally.

Chosen over STRING for this MVP: Enrichr is a two-call REST API with no
auth and returns ready-to-rank enrichment terms directly, which pairs
naturally with marker gene lists the fixed core already computes (e.g.
`rank_genes_groups` per cluster). STRING's interaction-network response is
richer but needs more client-side interpretation to turn into a short,
LLM-ready summary -- more surface area for an MVP tool.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import requests

from multiome_agent.config import REPO_ROOT
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

CACHE_DIR = REPO_ROOT / "cache" / "enrichr"
ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"
DEFAULT_LIBRARY = "GO_Biological_Process_2023"
TOP_N = 10


def _cache_key(genes: list[str], library: str) -> str:
    payload = json.dumps({"genes": sorted(genes), "library": library})
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def enrich_gene_set(genes: list[str], library: str = DEFAULT_LIBRARY, top_n: int = TOP_N) -> list[dict]:
    """Return up to `top_n` enriched terms for `genes` against `library`.

    Each result: {term, p_value, adjusted_p_value, genes_overlap}.
    """
    cache_file = CACHE_DIR / f"{_cache_key(genes, library)}.json"
    if cache_file.exists():
        logger.info("Enrichr cache hit for %d genes / %s", len(genes), library)
        return json.loads(cache_file.read_text())["results"]

    logger.info("Enrichr cache miss, querying for %d genes / %s", len(genes), library)
    add_resp = requests.post(
        f"{ENRICHR_BASE}/addList",
        files={"list": (None, "\n".join(genes)), "description": (None, "multiome_agent query")},
        timeout=15,
    )
    add_resp.raise_for_status()
    user_list_id = add_resp.json()["userListId"]

    enrich_resp = requests.get(
        f"{ENRICHR_BASE}/enrich",
        params={"userListId": user_list_id, "backgroundType": library},
        timeout=15,
    )
    enrich_resp.raise_for_status()
    rows = enrich_resp.json().get(library, [])

    # Enrichr row format: [rank, term, p_value, z_score, combined_score,
    # overlapping_genes, adjusted_p_value, ...]
    results = [
        {
            "term": row[1],
            "p_value": row[2],
            "adjusted_p_value": row[6],
            "genes_overlap": row[5],
        }
        for row in sorted(rows, key=lambda r: r[6])[:top_n]
    ]

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps({"genes": genes, "library": library, "results": results}, indent=2))
    return results
