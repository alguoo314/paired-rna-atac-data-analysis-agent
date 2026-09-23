"""PubMed E-utilities literature search, cached locally.

Two-step retrieval, matching how literature RAG actually needs to work: 1)
`search_pubmed` (esearch + esummary) finds CANDIDATE papers by keyword --
titles only, cheap, broad. 2) `fetch_pubmed_abstracts` (efetch, real abstract
text) reads the actual content of a short, already-narrowed list of PMIDs.
A claim grounded only in a title match is "direct keyword searching"; a
claim grounded in what an abstract's real text says is retrieval-augmented.
Kept as two separate tools/calls rather than always fetching full abstracts
for every search hit -- that would multiply token cost by the width of the
initial keyword search for no benefit when most hits get discarded anyway.
"""

from __future__ import annotations

import hashlib
import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

from multiome_agent.config import REPO_ROOT
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

CACHE_DIR = REPO_ROOT / "cache" / "pubmed"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_PARAMS = {"tool": "multiome_agent", "email": "multiome-agent@example.com"}
TOP_N = 5
MAX_ABSTRACT_CHARS = 2000  # per-abstract cap, same "truncate results" spirit as TOP_N

# NCBI throttles unauthenticated E-utilities calls to ~3 req/sec and returns
# 429 above that -- caught for real while testing `fetch_pubmed_abstracts`
# with a cleared cache (3 rapid uncached calls, no delay, immediate 429).
# The new RAG-style checklist/novelty steps make several sequential calls
# per report run, so a bare `raise_for_status()` would make the whole report
# generator flaky. Simple fixed-backoff retry, not exponential -- NCBI's
# limit is a steady rate, not a load spike to back off hard from. Also
# retries on 5xx: caught for real during an Opus report run -- two
# consecutive real "500 Internal Server Error" responses from NCBI itself
# (not a rate limit), which a real agent handled gracefully as a tool error
# but wasted turns/cost on retries the model had to improvise itself.
_MAX_RETRIES = 4
_RETRY_DELAY_S = 1.0
_RETRYABLE_STATUS = {429, 500, 502, 503, 504}


def _get_with_retry(url: str, params: dict, timeout: int) -> requests.Response:
    for attempt in range(_MAX_RETRIES):
        resp = requests.get(url, params=params, timeout=timeout)
        if resp.status_code not in _RETRYABLE_STATUS:
            resp.raise_for_status()
            return resp
        logger.warning(
            "PubMed E-utilities returned %d (attempt %d/%d), retrying",
            resp.status_code, attempt + 1, _MAX_RETRIES,
        )
        time.sleep(_RETRY_DELAY_S)
    resp.raise_for_status()  # exhausted retries -- raise the real error rather than returning a stale response
    return resp


def _cache_key(query: str) -> str:
    return hashlib.sha256(query.encode()).hexdigest()[:16]


def _cache_key_ids(pmids: list[str]) -> str:
    return hashlib.sha256(",".join(sorted(pmids)).encode()).hexdigest()[:16]


def search_pubmed(query: str, top_n: int = TOP_N) -> list[dict]:
    """Return up to `top_n` PubMed hits for `query` as [{pmid, title, journal, year}]."""
    cache_file = CACHE_DIR / f"{_cache_key(query)}.json"
    if cache_file.exists():
        logger.info("PubMed cache hit for query=%r", query)
        return json.loads(cache_file.read_text())["results"]

    logger.info("PubMed cache miss, querying E-utilities for query=%r", query)
    search_resp = _get_with_retry(
        f"{EUTILS_BASE}/esearch.fcgi",
        params={**TOOL_PARAMS, "db": "pubmed", "term": query, "retmax": top_n, "retmode": "json"},
        timeout=15,
    )
    pmids = search_resp.json()["esearchresult"]["idlist"]

    results = []
    if pmids:
        summary_resp = _get_with_retry(
            f"{EUTILS_BASE}/esummary.fcgi",
            params={**TOOL_PARAMS, "db": "pubmed", "id": ",".join(pmids), "retmode": "json"},
            timeout=15,
        )
        summary = summary_resp.json()["result"]
        for pmid in pmids:
            doc = summary.get(pmid, {})
            results.append({
                "pmid": pmid,
                "title": doc.get("title", ""),
                "journal": doc.get("fulljournalname", doc.get("source", "")),
                "year": (doc.get("pubdate", "") or "")[:4],
            })

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps({"query": query, "results": results}, indent=2))
    return results


def fetch_pubmed_abstracts(pmids: list[str]) -> dict[str, dict]:
    """Real abstract text for a small, already-narrowed list of PMIDs (via
    efetch, XML mode -- parsed per-article rather than plain text so
    multiple abstracts in one response can't bleed into each other).
    Returns {pmid: {title, abstract, journal, year}}; a PMID with no
    abstract on record gets an explicit "(no abstract available)" rather
    than being silently dropped, so the agent knows to look elsewhere for
    that specific claim rather than assuming it just wasn't returned.
    """
    if not pmids:
        return {}
    cache_file = CACHE_DIR / f"abstracts_{_cache_key_ids(pmids)}.json"
    if cache_file.exists():
        logger.info("PubMed abstract cache hit for pmids=%r", pmids)
        return json.loads(cache_file.read_text())

    logger.info("PubMed abstract cache miss, fetching for pmids=%r", pmids)
    resp = _get_with_retry(
        f"{EUTILS_BASE}/efetch.fcgi",
        params={**TOOL_PARAMS, "db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "xml"},
        timeout=20,
    )
    root = ET.fromstring(resp.content)

    results: dict[str, dict] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        if pmid_el is None or not pmid_el.text:
            continue
        pmid = pmid_el.text
        title_el = article.find(".//ArticleTitle")
        title = "".join(title_el.itertext()).strip() if title_el is not None else ""
        abstract = " ".join("".join(el.itertext()) for el in article.findall(".//AbstractText")).strip()
        if len(abstract) > MAX_ABSTRACT_CHARS:
            abstract = abstract[:MAX_ABSTRACT_CHARS] + "... (truncated)"
        journal_el = article.find(".//Journal/Title")
        year_el = article.find(".//JournalIssue/PubDate/Year")
        results[pmid] = {
            "title": title,
            "abstract": abstract or "(no abstract available)",
            "journal": journal_el.text if journal_el is not None else "",
            "year": year_el.text if year_el is not None else "",
        }

    # PMIDs the response didn't include an article for at all (rare, but
    # possible for a bad/retracted ID) -- explicit rather than silently missing.
    for pmid in pmids:
        results.setdefault(pmid, {"title": "", "abstract": "(no record found for this PMID)", "journal": "", "year": ""})

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps(results, indent=2))
    return results
