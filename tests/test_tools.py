"""Tests for the literature (PubMed) and database (Enrichr) tools.

Each makes one real network call (free public APIs, no auth) to confirm
the wrappers actually work, plus a cache-hit test following the pattern in
test_loader.py's cache round-trip test.
"""

import pytest

from multiome_agent.tools import enrichr, pubmed


def test_search_pubmed_real_call_returns_results():
    results = pubmed.search_pubmed("SPI1 PU.1 monocyte differentiation")
    assert 1 <= len(results) <= 5
    assert all({"pmid", "title", "journal", "year"} <= r.keys() for r in results)
    assert all(r["pmid"].isdigit() for r in results)


def test_search_pubmed_cache_hit(tmp_path, monkeypatch):
    monkeypatch.setattr(pubmed, "CACHE_DIR", tmp_path)
    query = "test cache query xyz123"

    first = pubmed.search_pubmed(query)
    cache_files = list(tmp_path.glob("*.json"))
    assert len(cache_files) == 1

    def _boom(*a, **k):
        raise AssertionError("cache miss: requests.get called when it shouldn't have been")

    monkeypatch.setattr(pubmed.requests, "get", _boom)
    second = pubmed.search_pubmed(query)
    assert second == first


def test_enrich_gene_set_real_call_returns_results():
    results = enrichr.enrich_gene_set(["CD14", "LYZ", "FCN1", "S100A8", "S100A9"])
    assert 1 <= len(results) <= 10
    assert all({"term", "p_value", "adjusted_p_value", "genes_overlap"} <= r.keys() for r in results)
    assert results == sorted(results, key=lambda r: r["adjusted_p_value"])


def test_enrich_gene_set_cache_hit(tmp_path, monkeypatch):
    monkeypatch.setattr(enrichr, "CACHE_DIR", tmp_path)
    genes = ["CD3E", "CD3D", "CD247"]

    first = enrichr.enrich_gene_set(genes)
    cache_files = list(tmp_path.glob("*.json"))
    assert len(cache_files) == 1

    def _boom(*a, **k):
        raise AssertionError("cache miss: requests called when it shouldn't have been")

    monkeypatch.setattr(enrichr.requests, "get", _boom)
    monkeypatch.setattr(enrichr.requests, "post", _boom)
    second = enrichr.enrich_gene_set(genes)
    assert second == first
