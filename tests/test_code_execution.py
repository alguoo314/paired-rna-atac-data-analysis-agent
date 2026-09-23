"""Tests for the local sandboxed code-execution tool (Stretch phase, Part B:
agent-written analysis code). Security-relevant: the denial tests matter as
much as the happy path.
"""

from __future__ import annotations

from multiome_agent.tools.code_execution import run_analysis_code

DATA = {"cluster_sizes": [100, 200, 150], "median_mito": [5.0, 9.9, 3.2]}


def test_happy_path_computes_real_correct_result():
    out = run_analysis_code('result = sum(data["cluster_sizes"])', DATA)
    assert "result = 450" in out


def test_happy_path_pandas_computes_real_correct_result():
    out = run_analysis_code('df = pd.DataFrame(data); result = float(df["median_mito"].mean())', DATA)
    assert "result = 6.033" in out


def test_print_output_is_captured():
    out = run_analysis_code('print("hello", 1 + 1)', DATA)
    assert "hello 2" in out


def test_rejects_import():
    out = run_analysis_code('import os; result = os.listdir(".")', DATA)
    assert out.startswith("Rejected")


def test_rejects_dunder_access():
    out = run_analysis_code("result = ().__class__.__bases__[0].__subclasses__()", DATA)
    assert out.startswith("Rejected")


def test_rejects_open():
    out = run_analysis_code('f = open("/etc/passwd"); result = f.read()', DATA)
    assert out.startswith("Rejected")


def test_rejects_exec_and_eval():
    assert run_analysis_code('exec("result = 1")', DATA).startswith("Rejected")
    assert run_analysis_code('result = eval("1+1")', DATA).startswith("Rejected")


def test_runtime_error_is_caught_not_raised():
    out = run_analysis_code("result = 1 / 0", DATA)
    assert "ZeroDivisionError" in out


def test_no_access_to_anything_outside_provided_namespace():
    # `namespace_data` (only `data`/`pd`/`np`) is the whole world -- nothing
    # else should be reachable, e.g. no access to this test's own DATA
    # object by any name other than what was explicitly passed as `data`.
    out = run_analysis_code("result = list(dict(__builtins__=1).keys())", {})
    assert out.startswith("Rejected")  # caught by the `__` denylist pattern


def test_empty_namespace_data_still_works():
    out = run_analysis_code("result = 1 + 1", {})
    assert "result = 2" in out


def test_rejects_pandas_numpy_file_io():
    # `pd`/`np` in the namespace are the real modules, not restricted APIs --
    # `open(` alone does NOT stop `pd.read_csv`/`np.load` from reading
    # arbitrary local files (verified directly against a real file before
    # this test was added: pd.read_csv("/etc/passwd", ...) succeeded and
    # returned real file contents). This is the single most important test
    # in this file given the project's core privacy design (no raw/local
    # data should ever be reachable from agent-written code) -- a gap here
    # would let the agent read the shareseq-multi-cell-lines dataset's own files directly.
    assert run_analysis_code('result = pd.read_csv("/etc/passwd").shape', {}).startswith("Rejected")
    assert run_analysis_code('result = np.load("/etc/passwd")', {}).startswith("Rejected")
    assert run_analysis_code('pd.DataFrame({"a":[1]}).to_csv("/tmp/leak.csv")', {}).startswith("Rejected")
