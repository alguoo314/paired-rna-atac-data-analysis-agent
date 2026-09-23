"""Automated grounding checks on an agent run's final answer: does every
number and every cited PMID actually trace back to a real tool result from
that same run? CLAUDE.md: "Every number in the report matches computed
tables. Every cited PMID / database ID exists."

Scope: this checks internal consistency (the model isn't citing a number or
PMID it never actually received from a tool this run) -- it does NOT
separately re-verify a PMID exists on PubMed itself, since a PMID that came
back from a real `search_pubmed` call in this same run is, by construction,
real (E-utilities wouldn't have returned it otherwise). A PMID appearing in
the answer that was never in any tool result is flagged as ungrounded
regardless of whether it happens to be a real PMID elsewhere -- the point is
that *this report* didn't earn the right to cite it.
"""

from __future__ import annotations

import math
import re

# Standalone decimal/scientific-notation numbers, e.g. 0.595, 3.5e-49.
# Requires a digit before or after the decimal point so bare words like
# "2026" alone aren't excluded, but a lone year is still a "number" by this
# regex; the grounding check only flags numbers that DON'T appear anywhere
# in the tool results, so an incidental year mentioned in prose that also
# happens to appear in some tool result text won't false-positive, and one
# that doesn't will -- authors should avoid bare incidental numbers in the
# answer, which is a reasonable thing to ask of a grounded scientific
# report anyway.
_NUMBER_RE = re.compile(r"-?\d+\.\d+(?:[eE][+-]?\d+)?|-?\d+[eE][+-]?\d+")

# LLM answers commonly render scientific notation as Unicode ("3.5 × 10⁻⁴⁹")
# rather than ASCII ("3.5e-49") -- without special-casing this, the plain
# regex above matches just the "3.5" mantissa fragment, which then fails to
# ground-check against a tool result like "3.5279961399705646e-49" (rounding
# 3.53e-49 to 1 decimal gives 0.0, not 3.5) and gets wrongly flagged as
# ungrounded. Match the whole Unicode-notation span first so it's evaluated
# as one number (3.5e-49), not decomposed into a misleading fragment.
_SUPERSCRIPT_TRANS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
_SCI_UNICODE_RE = re.compile(r"-?\d+(?:\.\d+)?\s*[×x]\s*10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)")

_PMID_RE = re.compile(r"PMID[:\s]*(\d{5,9})", re.IGNORECASE)

# Two distinct JSON shapes carry a real PMID in tool output: search_pubmed's
# list-of-records (`{"pmid": "123", ...}`) and fetch_pubmed_abstracts' dict
# keyed BY pmid (`{"123": {"title": ..., "abstract": ...}}` -- no literal
# "pmid" field at all). Missing the second shape would flag every PMID the
# agent only ever saw via fetch_pubmed_abstracts as ungrounded, even though
# it's a real tool-sourced PMID -- checked directly against that function's
# actual return shape, not assumed.
_TOOL_PMID_FIELD_RE = re.compile(r'"pmid":\s*"(\d+)"')
_TOOL_PMID_KEY_RE = re.compile(r'"(\d{5,9})":\s*\{')


def _extract_numbers(text: str) -> list[tuple[str, float]]:
    """Returns (display_string, parsed_float) pairs, preferring a full
    Unicode-scientific-notation match over the plain-number match it would
    otherwise be decomposed into."""
    found: list[tuple[int, int, str, float]] = []
    for m in _SCI_UNICODE_RE.finditer(text):
        exponent = int(m.group(1).translate(_SUPERSCRIPT_TRANS))
        mantissa = float(m.group(0).split("×" if "×" in m.group(0) else "x")[0].strip())
        found.append((m.start(), m.end(), m.group(0), mantissa * (10 ** exponent)))

    covered = [(s, e) for s, e, _, _ in found]
    for m in _NUMBER_RE.finditer(text):
        if any(s <= m.start() < e for s, e in covered):
            continue  # already part of a Unicode-sci-notation match
        found.append((m.start(), m.end(), m.group(0), float(m.group(0))))

    found.sort(key=lambda t: t[0])
    return [(disp, val) for _, _, disp, val in found]


_REL_TOL = 0.01  # 1% -- generous enough for typical display-rounding of a
# tool's raw float (e.g. citing 0.59 for a tool's 0.5948659415096681, ~0.8%
# off), tight enough to catch a genuinely different value (a value 2%+ off
# is a different number, not a rounding choice). Relative (not decimal-place
# exact) comparison because it handles Unicode scientific notation and plain
# decimals uniformly, without needing to define "how many decimal places"
# means for a citation like "3.5 x 10^-49".


def _cited_number_is_grounded(cited_f: float, tool_numbers: list[float]) -> bool:
    return any(
        (cited_f == 0 and t == 0) or math.isclose(cited_f, t, rel_tol=_REL_TOL)
        for t in tool_numbers
    )


def check_grounding(result) -> dict:
    """`result` is an `AgentRunResult` (or anything with `.answer` and
    `.tool_calls`, each tool_call a dict with a `result_summary` string).
    Returns {numbers: [...], pmids: [...], all_grounded: bool} where each
    entry records the cited value and whether it was found in this run's
    tool results.
    """
    # Prefer the full untruncated "result" field (added once tool outputs
    # started regularly exceeding the display-truncated "result_summary",
    # e.g. multi-abstract fetch_pubmed_abstracts calls) -- fall back to
    # "result_summary" for any older cached eval result that predates it.
    tool_text = " ".join(str(tc.get("result", tc.get("result_summary", ""))) for tc in result.tool_calls)
    tool_numbers = [val for _, val in _extract_numbers(tool_text)]

    number_findings = []
    for cited, cited_f in _extract_numbers(result.answer):
        grounded = _cited_number_is_grounded(cited_f, tool_numbers)
        number_findings.append({"cited": cited, "grounded": grounded})

    pmid_findings = []
    cited_pmids = set(_PMID_RE.findall(result.answer))
    tool_pmids = set(_TOOL_PMID_FIELD_RE.findall(tool_text)) | set(_TOOL_PMID_KEY_RE.findall(tool_text))
    for pmid in cited_pmids:
        pmid_findings.append({"pmid": pmid, "grounded": pmid in tool_pmids})

    all_grounded = all(n["grounded"] for n in number_findings) and all(p["grounded"] for p in pmid_findings)
    return {"numbers": number_findings, "pmids": pmid_findings, "all_grounded": all_grounded}
