"""Ties the fault injector to a real agent quality judgment, scored against
ground truth. Originally built for a per-model Haiku eval + a 4-model
comparison table (Day 2 step 7); that comparison table is retired (see
CLAUDE.md's Evaluation §8) in favor of the comprehensive report
(`agent/report_generator.py`), which reuses `build_scenarios`/`_run_scenario`
directly, one run per scenario, no consistency repeats. `run_model_eval`
(with its optional consistency-repeat wrapper) is kept as a lighter-weight
single-model entry point.

Design choice: every scenario's `run_agent` call gets the CLEAN
`agent_fixed_core_mdata` as `mdata` (so `tf_motif_correlation`/`enrich_gene_set`
keep working sensibly -- none of the 3 faults modify the peak count matrix or
chromVAR deviations those tools read from) but a scenario-specific
`qc_summary` composed from step 5's already-validated `score_*` deltas. Only
`get_qc_summary`'s output differs per scenario; the agent has to actually
call that tool and read real numbers to notice anything, not have numbers
smuggled into the system prompt where they'd escape the grounding check.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from anndata import AnnData
from mudata import MuData

from multiome_agent.agent.loop import run_agent
from multiome_agent.agent.qc_summary import fixed_core_summary, format_qc_summary
from multiome_agent.config import REPO_ROOT
from multiome_agent.core.rna_qc import run_rna_qc
from multiome_agent.data.loader import RAW_DIR
from multiome_agent.eval.fault_injection import (
    clean_control,
    downsample_atac_fragments,
    inject_doublets,
    shuffle_rna_atac_pairing,
)
from multiome_agent.eval.grounding import check_grounding
from multiome_agent.eval.scoring import score_doublets, score_downsampling, score_shuffled_pairing
from multiome_agent.logging_utils import get_logger

logger = get_logger(__name__)

FRAGMENTS_FILE = RAW_DIR / "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
DOWNSAMPLE_SCRATCH_DIR = REPO_ROOT / "data" / "processed" / "eval_scratch"

QUESTION = (
    "Assess the quality of this processed multiome dataset using the QC summary "
    "tool. Does anything look wrong? If so, describe the specific problem and how "
    "confident you are; if not, say the data looks clean."
)

# Heuristic free-text classifiers, not exact ground-truth matches -- a
# genuinely open question ("does this look wrong") doesn't have a single
# correct string. Documented as a heuristic in PROGRESS.md, not silently
# treated as exact.
PROBLEM_WORDS = [
    "wrong", "problem", "issue", "concern", "unusual", "abnormal", "unexpected",
    "red flag", "elevated", "degraded", "suspicious", "poor quality", "flagged",
    "inconsistent", "discrepanc", "mismatch",
]
# Cues that negate a problem-word hit immediately before it, e.g. "no major
# red flags" or "not a problem" -- caught by an actual test failing on real
# clean-summary phrasing ("This dataset looks reasonably healthy overall. No
# major red flags."), not written defensively up front. A ~20-char lookback
# window is a blunt instrument (won't catch negation further back in the
# sentence), but it's enough for the common "no/not X" pattern LLM answers
# actually use, and this is documented as a heuristic classifier, not an
# exact one.
_NEGATION_CUES = ("no ", "not ", "n't ", "without ", "none ", "nothing")
_NEGATION_WINDOW = 20


def _classify_detected(answer_text: str) -> bool:
    text = answer_text.lower()
    for word in PROBLEM_WORDS:
        start = 0
        while (idx := text.find(word, start)) != -1:
            window = text[max(0, idx - _NEGATION_WINDOW):idx]
            if not any(cue in window for cue in _NEGATION_CUES):
                return True
            start = idx + len(word)
    return False


FAULT_KEYWORDS = {
    "shuffled_rna_atac_pairing": ["pairing", "mismatch", "cluster agreement", "cross-modal", " ari", "shuffl", "misalign", "correspond"],
    "atac_downsampling": ["fragment", "depth", "downsampl", "sequencing depth", "coverage", "tss enrichment", "shallow"],
    "injected_doublets": ["doublet", "multiplet"],
}


def _classify_diagnosis_matches(answer_text: str, fault_type: str) -> bool | None:
    if fault_type == "clean_control":
        return None
    return any(kw in answer_text.lower() for kw in FAULT_KEYWORDS[fault_type])


@dataclass
class ScenarioResult:
    fault_type: str
    severity: float
    description: str
    answer: str
    detected: bool
    diagnosis_matches: bool | None
    false_alarm: bool
    cost: float
    turns: int
    grounding: dict
    reasoning_trail_path: str
    repeat_index: int = 0
    refused: bool = False
    refusal_detail: str | None = None


def _qc_summary_clean(clean_mdata) -> str:
    return format_qc_summary(fixed_core_summary(clean_mdata))


def _qc_summary_for_scenario(clean_mdata, fault_type: str, fault_score: dict | None) -> str:
    """Builds ONE coherent QC summary reflecting the (possibly faulted)
    dataset state -- never a "baseline vs. current" comparison. An earlier
    version stated the clean baseline alongside the faulted value directly
    in the tool result (e.g. "ARI is 0.140 (this dataset type's baseline ARI
    is 0.467)"), which handed the model the answer to "did anything change"
    for free rather than testing whether it can judge a single QC snapshot
    against its own knowledge of typical ranges -- caught by the project
    owner reading the committed report, not found independently. Fixed by
    overriding only the fault-affected fields in the summary dict and
    rendering through the SAME `format_qc_summary` used for the clean case,
    so the faulted scenarios read as an ordinary QC summary a real analyst
    would receive, with no second (clean) number anywhere in the text.

    `injected_doublets` originally described the synthetic cells as "extra"
    alongside the "original" ones -- caught as the SAME category of leak by
    the project owner (a real doublet-contaminated dataset doesn't come with
    its droplets pre-labeled "real" vs. "extra" either): fixed to report ONE
    blended aggregate over the full augmented population, with no subgroup
    breakdown at all. The ATAC-side fields still reflect only the clean
    500-cell baseline (recomputing ATAC clustering/motifs for synthetic
    cells that have no real underlying fragments is out of scope here) --
    a known, named simplification, not silently glossed over.
    """
    summary = fixed_core_summary(clean_mdata)
    if fault_type == "clean_control":
        pass
    elif fault_type == "shuffled_rna_atac_pairing":
        summary = {**summary, "ari": fault_score["ari_faulted"], "spi1_motif_rho": fault_score["tf_motif_rho_faulted"]}
    elif fault_type == "atac_downsampling":
        summary = {
            **summary,
            "median_atac_fragments": fault_score["n_fragment_median_faulted"],
            "median_tss_enrichment": fault_score["tss_enrichment_median_faulted"],
            "median_nucleosome_signal": fault_score["nucleosome_signal_median_faulted"],
        }
    elif fault_type == "injected_doublets":
        summary = {
            **summary,
            "n_cells": fault_score["blended_n_cells"],
            "median_genes_per_cell": fault_score["blended_median_genes_per_cell"],
            "median_umis_per_cell": fault_score["blended_median_umis_per_cell"],
            "median_pct_mt": fault_score["blended_median_pct_mt"],
            "predicted_doublet_rate": fault_score["blended_predicted_doublet_rate"],
            "median_doublet_score": fault_score["blended_median_doublet_score"],
        }
    else:
        raise ValueError(f"unknown fault_type {fault_type!r}")
    return format_qc_summary(summary)


def build_scenarios(clean_mdata) -> list[dict]:
    """Returns a list of {fault_type, severity, description, qc_summary}."""
    scenarios = []

    _, cc_meta = clean_control(clean_mdata)
    scenarios.append({
        "fault_type": "clean_control", "severity": 0.0, "description": cc_meta.description,
        "qc_summary": _qc_summary_for_scenario(clean_mdata, "clean_control", None),
    })

    shuffled_mdata, sh_meta = shuffle_rna_atac_pairing(clean_mdata, fraction=0.5, seed=1)
    sh_score = score_shuffled_pairing(clean_mdata, shuffled_mdata)
    scenarios.append({
        "fault_type": "shuffled_rna_atac_pairing", "severity": 0.5, "description": sh_meta.description,
        "qc_summary": _qc_summary_for_scenario(clean_mdata, "shuffled_rna_atac_pairing", sh_score),
    })

    atac_barcodes = list(clean_mdata["atac"].obs_names)
    DOWNSAMPLE_SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    ds_path, ds_meta = downsample_atac_fragments(
        FRAGMENTS_FILE, atac_barcodes, keep_frac=0.2, out_dir=DOWNSAMPLE_SCRATCH_DIR, seed=0
    )
    logger.info("Running real snapatac2 rerun for the downsampling scenario (~17 min)...")
    ds_score = score_downsampling(clean_mdata["atac"], atac_barcodes, FRAGMENTS_FILE, ds_path, keep_frac=0.2)
    scenarios.append({
        "fault_type": "atac_downsampling", "severity": 0.8, "description": ds_meta.description,
        "qc_summary": _qc_summary_for_scenario(clean_mdata, "atac_downsampling", ds_score),
    })

    # `inject_doublets` sums `.X` directly and explicitly requires RAW counts
    # (its own docstring: summing log-normalized/TF-IDF values isn't a valid
    # doublet simulation) -- but `clean_mdata` here is the agent's fixed-core
    # cache, already through `run_fixed_core`'s clustering step, which
    # overwrites `.X` with normalized/TF-IDF values and stashes the raw
    # counts in `.layers["counts"]` instead (see clustering.py). Build a
    # raw-counts-only MuData for this specific call rather than passing
    # `clean_mdata` directly -- caught by checking clustering.py's actual
    # behavior before running the expensive real eval, not by assuming
    # `clean_mdata.X` was still raw.
    raw_rna = AnnData(
        X=clean_mdata["rna"].layers["counts"].copy(), obs=clean_mdata["rna"].obs.copy(), var=clean_mdata["rna"].var.copy()
    )
    raw_atac = AnnData(
        X=clean_mdata["atac"].layers["counts"].copy(), obs=clean_mdata["atac"].obs.copy(), var=clean_mdata["atac"].var.copy()
    )
    raw_mdata = MuData({"rna": raw_rna, "atac": raw_atac})
    doublet_mdata, db_meta = inject_doublets(raw_mdata, doublet_rate=0.1, seed=0)
    rna_faulted = doublet_mdata["rna"].copy()
    run_rna_qc(rna_faulted)
    db_score = score_doublets(rna_faulted.obs, db_meta.ground_truth["doublet_obs_names"])
    # Blended aggregates over the FULL augmented population (real + synthetic
    # together), not a "real vs. extra" breakdown -- naming a subgroup as
    # "extra"/"core" hands the model the same kind of answer-shaped hint as
    # the baseline-leak bug fixed above (a real synthetic doublet is not
    # labeled as such in a real dataset either). `rna_faulted.obs` already
    # has genes/UMIs/mito/doublet QC computed over all 556 cells together.
    db_score = {
        **db_score,
        "blended_n_cells": rna_faulted.n_obs,
        "blended_median_genes_per_cell": float(rna_faulted.obs["n_genes_by_counts"].median()),
        "blended_median_umis_per_cell": float(rna_faulted.obs["total_counts"].median()),
        "blended_median_pct_mt": float(rna_faulted.obs["pct_counts_mt"].median()),
        "blended_predicted_doublet_rate": float(rna_faulted.obs["predicted_doublet"].mean()),
        "blended_median_doublet_score": float(rna_faulted.obs["doublet_score"].median()),
    }
    scenarios.append({
        "fault_type": "injected_doublets", "severity": 0.1, "description": db_meta.description,
        "qc_summary": _qc_summary_for_scenario(clean_mdata, "injected_doublets", db_score),
    })

    return scenarios


def _run_scenario(clean_mdata, scenario: dict, repeat_index: int = 0, model: str | None = None) -> ScenarioResult:
    result = run_agent(QUESTION, model=model, mdata=clean_mdata, qc_summary=scenario["qc_summary"])
    grounding = check_grounding(result)
    detected = _classify_detected(result.answer)
    is_clean = scenario["fault_type"] == "clean_control"
    return ScenarioResult(
        fault_type=scenario["fault_type"], severity=scenario["severity"], description=scenario["description"],
        answer=result.answer, detected=detected,
        diagnosis_matches=_classify_diagnosis_matches(result.answer, scenario["fault_type"]),
        false_alarm=(is_clean and detected), cost=result.estimated_cost_usd, turns=result.turn_count,
        grounding=grounding, reasoning_trail_path=result.reasoning_trail_path, repeat_index=repeat_index,
        refused=result.refused, refusal_detail=result.refusal_detail,
    )


def run_model_eval(model: str, scenarios: list[dict], clean_mdata, n_consistency_repeats: int = 0) -> dict:
    """Run the (already-built, reused) scenarios against any model. Building
    `scenarios` is expensive (one real ~17-min snapatac2 rerun for the
    downsampling scenario) -- callers comparing multiple models should build
    `scenarios` once via `build_scenarios` and pass the same list to every
    `run_model_eval` call, not rebuild per model.
    """
    scenario_results = [_run_scenario(clean_mdata, s, model=model) for s in scenarios]

    clean_scenario = scenarios[0]
    consistency_results = [
        _run_scenario(clean_mdata, clean_scenario, repeat_index=i + 1, model=model)
        for i in range(n_consistency_repeats)
    ]

    total_cost = sum(r.cost for r in scenario_results) + sum(r.cost for r in consistency_results)
    return {
        "model": model,
        "scenario_results": [asdict(r) for r in scenario_results],
        "consistency_results": [asdict(r) for r in consistency_results],
        "total_cost_usd": total_cost,
    }
