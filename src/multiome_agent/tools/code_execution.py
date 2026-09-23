"""Local, restricted Python execution over pre-aggregated summary data only.

Why not Anthropic's server-side `code_execution` tool: that requires
uploading data to Anthropic's own sandbox to be useful, which conflicts with
CLAUDE.md's "summaries, not raw data, go to the LLM" principle -- this is
especially important for the shareseq-multi-cell-lines dataset, where no per-cell data should
ever leave the local machine. This tool runs locally instead, and the ONLY
data ever placed in its namespace is whatever small aggregate dict the
caller explicitly builds ahead of time (e.g. per-cluster QC medians, a
handful of summary statistics) -- never a raw AnnData/MuData or per-cell
matrix, and callers must anonymize any identity strings before including
them, same convention as `eval/shareseq_fault_injection.py`.

This is a restricted namespace plus a string denylist, not a real sandbox:
it stops the obvious escape patterns (`import`, dunder access, file/network/
process APIs) but is not a formal security boundary against a genuinely
adversarial payload. Scoped for a cooperative LLM assistant investigating
its own already-aggregated summary data, not for running untrusted code.
"""

from __future__ import annotations

import builtins
import contextlib
import io

import numpy as np
import pandas as pd

_DENYLIST = [
    "import", "exec(", "eval(", "open(", "__", "os.", "sys.", "subprocess",
    "globals(", "locals(", "compile(", "getattr(", "setattr(",
    # `pd`/`np` are the real modules, not restricted APIs -- `open(` alone
    # doesn't stop `pd.read_csv(...)`/`np.load(...)` from reading arbitrary
    # local files (verified directly: pd.read_csv("/etc/passwd", ...)
    # succeeds without this). Denylist their I/O surface explicitly, both
    # read (exfiltrating local files, including the shareseq-multi-cell-lines dataset itself,
    # into the agent's context) and write (writing data out anywhere).
    "read_csv", "read_json", "read_excel", "read_parquet", "read_pickle",
    "read_hdf", "read_sql", "read_fwf", "read_table", "read_html",
    "read_xml", "read_orc", "read_feather", "read_stata", "read_sas",
    "read_spss", "read_gbq", "read_clipboard", "HDFStore", "ExcelFile",
    "ExcelWriter", "to_csv(", "to_pickle(", "to_excel(", "to_hdf(",
    "to_parquet(", "to_json(", "to_sql(", "to_clipboard(",
    "np.load(", "np.loadtxt(", "np.fromfile(", "np.genfromtxt(",
    "np.memmap(", "np.fromregex(", "np.save(", "np.savetxt(", "np.savez(",
]

_MAX_RESULT_CHARS = 2000

# Reference the real `builtins` module explicitly rather than the module-local
# `__builtins__` name, whose type (module vs. dict) differs between `__main__`
# and imported modules -- unambiguous either way.
_SAFE_BUILTINS = {
    name: getattr(builtins, name)
    for name in (
        "len", "range", "sum", "min", "max", "round", "abs", "sorted",
        "list", "dict", "set", "tuple", "float", "int", "str", "bool",
        "print", "enumerate", "zip", "map", "filter", "any", "all",
    )
}


def run_analysis_code(code: str, namespace_data: dict) -> str:
    """Executes `code` with only `pd`, `np`, and `data` (=`namespace_data`,
    a plain dict of pre-aggregated statistics) in scope. Returns captured
    stdout plus a `result = ...` line if the code assigned a `result`
    variable, or an error string -- never raises to the caller.
    """
    for pattern in _DENYLIST:
        if pattern in code:
            return f"Rejected: code contains disallowed pattern {pattern!r}."

    local_ns: dict = {"pd": pd, "np": np, "data": namespace_data}
    stdout = io.StringIO()
    try:
        with contextlib.redirect_stdout(stdout):
            exec(compile(code, "<agent_analysis_code>", "exec"), {"__builtins__": _SAFE_BUILTINS}, local_ns)
    except Exception as e:
        return f"Execution error: {type(e).__name__}: {e}"

    out = stdout.getvalue()
    if "result" in local_ns:
        out += f"\nresult = {local_ns['result']!r}"
    if len(out) > _MAX_RESULT_CHARS:
        out = out[:_MAX_RESULT_CHARS] + "... (truncated)"
    return out or "(no output -- assign to a variable named `result` or use print())"
