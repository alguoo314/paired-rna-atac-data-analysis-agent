"""Tests for the agent loop: max-turns guard (mocked, no real API calls),
reasoning-trail format, and one real end-to-end run on Haiku.
"""

import os

import pytest

from multiome_agent.agent import loop as loop_mod


class _FakeTextBlock:
    type = "text"

    def __init__(self, text):
        self.text = text


class _FakeToolUseBlock:
    type = "tool_use"

    def __init__(self, name, tool_input, block_id):
        self.name = name
        self.input = tool_input
        self.id = block_id


class _FakeUsage:
    input_tokens = 10
    output_tokens = 10
    cache_read_input_tokens = 0
    cache_creation_input_tokens = 0


class _FakeResponse:
    def __init__(self, content, stop_reason):
        self.content = content
        self.stop_reason = stop_reason
        self.usage = _FakeUsage()


class _AlwaysToolUseClient:
    """Stub client that always asks to call a nonexistent tool -- never
    reaches end_turn, so the loop must rely on the max-turns guard alone.
    Nonexistent tool name means `_execute_tool` never touches `mdata`.
    """

    call_count = 0

    def __init__(self, api_key=None):
        self.messages = self

    def create(self, **kwargs):
        _AlwaysToolUseClient.call_count += 1
        return _FakeResponse(
            content=[
                _FakeTextBlock(f"thinking on call {_AlwaysToolUseClient.call_count}"),
                _FakeToolUseBlock("nonexistent_tool", {}, f"id{_AlwaysToolUseClient.call_count}"),
            ],
            stop_reason="tool_use",
        )


def test_max_turns_guard_stops_exactly_at_limit(monkeypatch, tmp_path):
    _AlwaysToolUseClient.call_count = 0
    monkeypatch.setattr(loop_mod.anthropic, "Anthropic", _AlwaysToolUseClient)
    monkeypatch.setattr(loop_mod, "get_agent_fixed_core", lambda: None)
    monkeypatch.setattr(loop_mod, "RUNS_DIR", tmp_path)
    monkeypatch.setattr(loop_mod, "get_anthropic_api_key", lambda: "fake-key-not-used")

    result = loop_mod.run_agent("infinite question", model="claude-haiku-4-5", max_turns=3)

    assert result.hit_max_turns is True
    assert result.turn_count == 3
    assert _AlwaysToolUseClient.call_count == 3
    assert "max-turns" in result.answer.lower()


def test_reasoning_trail_is_readable_markdown(tmp_path, monkeypatch):
    monkeypatch.setattr(loop_mod, "RUNS_DIR", tmp_path)
    trail = loop_mod.ReasoningTrail("what is SPI1?", "claude-haiku-4-5")
    trail.log_turn(
        1, ["I should check the TF-motif correlation."],
        [{"name": "tf_motif_correlation", "input": {"gene": "SPI1"}, "result_summary": '{"gene": "SPI1", "spearman_rho": 0.6}'}],
    )
    trail.log_final("SPI1 expression correlates with its motif deviation.", False, 0.0012, 1)

    content = trail.path.read_text()
    assert "what is SPI1?" in content
    assert "tf_motif_correlation" in content
    assert "SPI1 expression correlates" in content
    assert "Turns:** 1" in content
    assert trail.path.suffix == ".md"


class _EndsWithoutRequiredToolClient:
    """Stub that ends on `end_turn` without calling `require_tool_call` on its
    first response, then (once nudged with a forced `tool_choice`) calls it
    on the second. Reproduces the real TCF7L2-judging run: the model
    spiralled through literature search and ended with prose instead of
    calling `record_judger_verdict`.
    """

    def __init__(self, api_key=None):
        self.messages = self
        self.calls = []

    def create(self, **kwargs):
        self.calls.append(kwargs)
        if len(self.calls) == 1:
            return _FakeResponse(content=[_FakeTextBlock("my conclusion, no tool call")], stop_reason="end_turn")
        if len(self.calls) == 2:
            assert kwargs.get("tool_choice") == {"type": "tool", "name": "record_judger_verdict"}
            return _FakeResponse(
                content=[_FakeToolUseBlock("record_judger_verdict", {"verdict": "struck_down"}, "id2")],
                stop_reason="tool_use",
            )
        # Third call: the required tool was already satisfied on call 2, so
        # the loop must accept this end_turn instead of nudging again.
        return _FakeResponse(content=[_FakeTextBlock("done")], stop_reason="end_turn")


def test_require_tool_call_nudges_once_then_succeeds(monkeypatch, tmp_path):
    client = _EndsWithoutRequiredToolClient()
    monkeypatch.setattr(loop_mod.anthropic, "Anthropic", lambda api_key=None: client)
    monkeypatch.setattr(loop_mod, "get_agent_fixed_core", lambda: None)
    monkeypatch.setattr(loop_mod, "RUNS_DIR", tmp_path)
    monkeypatch.setattr(loop_mod, "get_anthropic_api_key", lambda: "fake-key-not-used")
    monkeypatch.setattr(
        loop_mod, "_execute_tool",
        lambda name, inp, mdata, qc_summary, ns: (loop_mod.json.dumps(inp), False),
    )

    result = loop_mod.run_agent(
        "judge this finding", model="claude-haiku-4-5", max_turns=10,
        require_tool_call="record_judger_verdict",
    )

    assert len(client.calls) == 3
    assert any(tc["name"] == "record_judger_verdict" for tc in result.tool_calls)
    assert result.hit_max_turns is False


class _NeverCallsRequiredToolClient:
    """Even after being nudged with a forced `tool_choice`, this stub keeps
    ending on `end_turn` -- the nudge must fire at most once, not loop
    forever chasing a tool call that never comes.
    """

    def __init__(self, api_key=None):
        self.messages = self
        self.calls = []

    def create(self, **kwargs):
        self.calls.append(kwargs)
        return _FakeResponse(content=[_FakeTextBlock("still no tool call")], stop_reason="end_turn")


def test_require_tool_call_gives_up_after_one_nudge(monkeypatch, tmp_path):
    client = _NeverCallsRequiredToolClient()
    monkeypatch.setattr(loop_mod.anthropic, "Anthropic", lambda api_key=None: client)
    monkeypatch.setattr(loop_mod, "get_agent_fixed_core", lambda: None)
    monkeypatch.setattr(loop_mod, "RUNS_DIR", tmp_path)
    monkeypatch.setattr(loop_mod, "get_anthropic_api_key", lambda: "fake-key-not-used")

    result = loop_mod.run_agent(
        "judge this finding", model="claude-haiku-4-5", max_turns=10,
        require_tool_call="record_judger_verdict",
    )

    assert len(client.calls) == 2  # original attempt + exactly one nudge
    assert result.tool_calls == []
    assert result.hit_max_turns is False


@pytest.mark.skipif(not os.environ.get("ANTHROPIC_API_KEY"), reason="no ANTHROPIC_API_KEY available")
def test_real_end_to_end_run_on_haiku():
    result = loop_mod.run_agent(
        "In this PBMC dataset, does SPI1 expression correlate with its own motif's "
        "accessibility, and is that consistent with SPI1's known biology?",
        model="claude-haiku-4-5",
    )
    assert result.answer.strip() != ""
    assert result.turn_count >= 1
    assert result.estimated_cost_usd > 0
    assert os.path.exists(result.reasoning_trail_path)
    assert len(result.tool_calls) >= 1
