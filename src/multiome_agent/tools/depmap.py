"""Resolves a Broad Institute DepMap Model ID (format "ACH-XXXXXX") to its
real cell-line name and disease/lineage annotation, via EBI's public
Cellosaurus cell-line database (which cross-references DepMap IDs) --
DepMap's own portal (depmap.org) sits behind an interactive Cloudflare
bot-check that a plain HTTP client can't pass, so this project uses
Cellosaurus's REST API instead, verified directly against real ACH IDs
(ACH-000001 -> OVCAR-3, ACH-000004 -> HEL, ACH-000958 -> SW48) before
relying on it.

Deliberately no hardcoded ID->name table: this is a live lookup against a
real reference database, not a static mapping baked into this project's
code -- the agent is expected to notice an ACH-formatted value itself (e.g.
in a column literally named "Depmap") and decide to call this tool, not be
told to.
"""

from __future__ import annotations

import hashlib
import json

import requests

from multiome_agent.config import REPO_ROOT
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

CACHE_DIR = REPO_ROOT / "cache" / "depmap"
CELLOSAURUS_SEARCH_URL = "https://api.cellosaurus.org/search/cell-line"


def _cache_key(depmap_id: str) -> str:
    return hashlib.sha256(depmap_id.encode()).hexdigest()[:16]


def resolve_depmap_id(depmap_id: str) -> dict:
    """Returns {"found": bool, "cell_line_name": str|None, "disease": str|None,
    "species": str|None, "cellosaurus_accession": str|None}."""
    cache_file = CACHE_DIR / f"{_cache_key(depmap_id)}.json"
    if cache_file.exists():
        logger.info("DepMap ID cache hit for %r", depmap_id)
        return json.loads(cache_file.read_text())

    logger.info("DepMap ID cache miss, querying Cellosaurus for %r", depmap_id)
    resp = requests.get(
        CELLOSAURUS_SEARCH_URL,
        params={"q": f"dr:{depmap_id}", "fields": "id,name,ac,di,ox", "format": "json"},
        timeout=15,
    )
    resp.raise_for_status()
    hits = resp.json().get("Cellosaurus", {}).get("cell-line-list", [])

    if not hits:
        result = {"found": False, "cell_line_name": None, "disease": None, "species": None, "cellosaurus_accession": None}
    else:
        hit = hits[0]
        name = next((n["value"] for n in hit.get("name-list", []) if n.get("type") == "identifier"), None)
        disease = hit.get("disease-list", [{}])[0].get("label")
        species = hit.get("species-list", [{}])[0].get("label")
        accession = next((a["value"] for a in hit.get("accession-list", []) if a.get("type") == "primary"), None)
        result = {
            "found": True, "cell_line_name": name, "disease": disease,
            "species": species, "cellosaurus_accession": accession,
        }

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps(result, indent=2))
    return result
