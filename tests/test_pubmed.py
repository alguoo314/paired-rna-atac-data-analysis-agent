"""Tests for the two-step PubMed retrieval (search by keyword, then fetch
real abstract text for a narrowed PMID list) -- see tools/pubmed.py's
module docstring for why these are separate calls. Hits the real E-utilities
API (no API key required, matches this repo's other network-tool tests
which also aren't key-gated).
"""

from __future__ import annotations
from unittest.mock import MagicMock, patch

from multiome_agent.tools.pubmed import _get_with_retry, fetch_pubmed_abstracts, search_pubmed


def test_get_with_retry_retries_on_real_500_not_just_429():
    # Real bug: NCBI returned two consecutive genuine 500s (not rate-limit
    # 429s) during a real Opus report run -- the old code only retried 429
    # and would have raised immediately on a 500.
    ok_resp = MagicMock(status_code=200)
    err_resp = MagicMock(status_code=500)
    with patch("multiome_agent.tools.pubmed.requests.get", side_effect=[err_resp, err_resp, ok_resp]), \
         patch("multiome_agent.tools.pubmed.time.sleep"):
        result = _get_with_retry("http://x", {}, timeout=5)
    assert result is ok_resp
    ok_resp.raise_for_status.assert_called_once()


def test_search_pubmed_returns_titles_only():
    results = search_pubmed("SPI1 PU.1 myeloid chromatin accessibility")
    assert len(results) > 0
    for hit in results:
        assert hit["pmid"]
        assert hit["title"]
        assert "abstract" not in hit  # confirms this stays title-only, not the RAG step


def test_fetch_pubmed_abstracts_returns_real_content():
    hits = search_pubmed("SPI1 PU.1 myeloid chromatin accessibility")
    pmids = [h["pmid"] for h in hits[:2]]
    abstracts = fetch_pubmed_abstracts(pmids)
    assert set(abstracts.keys()) == set(pmids)
    for pmid in pmids:
        entry = abstracts[pmid]
        assert entry["title"]
        # Real abstract text, not just a restated title or empty placeholder.
        assert len(entry["abstract"]) > 100
        assert entry["abstract"] != "(no abstract available)"


def test_fetch_pubmed_abstracts_handles_missing_pmid_explicitly():
    result = fetch_pubmed_abstracts(["999999999999"])  # not a real PMID
    assert result["999999999999"]["abstract"] == "(no record found for this PMID)"


def test_fetch_pubmed_abstracts_empty_input():
    assert fetch_pubmed_abstracts([]) == {}
