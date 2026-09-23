"""Tests for DepMap Model ID resolution (tools/depmap.py) -- hits the real
Cellosaurus public API (no key required), verified against real, known
ACH IDs.
"""

from __future__ import annotations

from multiome_agent.tools.depmap import resolve_depmap_id


def test_resolve_known_depmap_id_returns_real_cell_line():
    result = resolve_depmap_id("ACH-000001")
    assert result["found"] is True
    assert result["cell_line_name"] == "OVCAR-3"
    assert "ovarian" in result["disease"].lower()


def test_resolve_unknown_depmap_id_returns_not_found():
    result = resolve_depmap_id("ACH-999999")
    assert result["found"] is False
    assert result["cell_line_name"] is None
