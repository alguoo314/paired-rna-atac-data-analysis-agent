# Multiome QC & Hypothesis Agent

An LLM agent that analyzes paired single-cell RNA + ATAC (multiome) data. It judges data
quality, determines what the data actually IS from evidence (never told directly —
cell type, cell line, or condition), builds a known-biology checklist grounded in
literature it actually reads (not just title-matches), and proposes candidate novel
findings that a separate, adversarial "Judger" agent tries to strike down. All numbers
(QC metrics, clustering, motif scores, correlations) are computed by validated Python
code; the LLM plans which analyses to run and interprets the results — it never invents
a number. A fault-injection benchmark measures whether the agent actually notices real,
planted data-quality problems, including a negative control checking it doesn't
hallucinate a relationship from data with no real cross-modal structure left.

This repo documents its own construction in detail in [`PROGRESS.md`](PROGRESS.md) —
every real number, every bug found and how, and every design decision, written as the
work happened. It's worth reading if you want to see what actually went wrong and how
it was diagnosed, not just the finished result.

## Architecture

```
agent-decided loader -> fixed core analyses -> agent loop (planner) -> 5-section report
  (inspects real file                              |
   structure; picks       +--------------+--------------+--------------+
   combined vs. per-      |              |              |              |
   modality strategy) analysis menu  literature RAG  identity/cross-  condition-axis
                      (validated       (search, then    modal tools     detection
                       wrapper)        fetch+read      (marker discovery,
                                       real abstracts)  cross_modal_marker_check)
```

**Fixed core** (always runs, deterministic, gives a stable baseline): RNA QC
(genes/UMIs/%mito/doublet scores via scrublet), ATAC QC (fragments/cell, TSS enrichment,
FRiP, nucleosome signal — computed directly from the real fragments file via
`snapatac2`/`muon`), Leiden clustering on both modalities + RNA-vs-ATAC cluster
agreement (ARI), gene activity scores, and real chromVAR-style motif deviation scores
(full JASPAR 2024 vertebrate motif set, GC/accessibility-matched background sampling,
actual genomic sequence from hg38).

**Loading is agent-mediated, not pre-specified**: a small dedicated agent call inspects
real file structure (shape, column names, feature-type values — never the matrix data
itself) and picks a named loading strategy via a schema-constrained tool call before
anything else runs. The actual parsing stays fully deterministic Python; only the
*choice* of which parser applies is agent-decided (`agent/loader_selector.py`).

**Agent loop**: a hand-rolled `client.messages.create` loop (not the SDK's beta tool
runner — a deliberate explainability choice, see PROGRESS.md Day 1 step 4) with a
max-turns guard and a human-readable reasoning-trail log per run
(`logs/agent_runs/*.md`, gitignored). Eleven tools: TF-expression-vs-own-motif-accessibility,
two-step PubMed RAG (search titles, then fetch and actually read real abstracts before
citing), Enrichr gene-set enrichment, identity-metadata detection (does this dataset's
own `.obs` already encode identity — by column name, or by value pattern, e.g. a Broad
Institute DepMap "ACH-XXXXXX" ID resolved via a live EBI Cellosaurus lookup) plus
marker-based fallback discovery on both the RNA and ATAC (gene-activity) side,
cross-modal cell-type-call validation, condition-axis detection + per-arm QC, and a
QC-summary tool.

**Identity is discovered, never told.** The system prompt does not say what tissue,
cell line, or condition a dataset is — that's exactly what the agent has to determine
from cluster markers and cross-modal validation before it can do anything
identity-specific (build a literature-grounded checklist, propose novel findings).
Verified on a real run with zero hints: the agent correctly reconstructed full PBMC
cell-type composition (monocyte/B/T/NK subsets, cross-modal-confirmed) from tools alone.

**Eval**: a fault injector corrupts copies of the clean dataset (shuffled RNA-ATAC
pairing, heavy ATAC downsampling, injected doublets, plus dataset-specific label/sample
faults — see the comprehensive report's §4), and a grounding checker verifies every
number and PMID an agent cites actually traces to a real tool result from that same run
(it has caught real fabricated citations, and a real motif-name-matching bug in this
project's own code — see PROGRESS.md).

## Results at a glance

The flagship deliverable is one comprehensive 5-section report per dataset — dataset
summary, gene activity/chromVAR + cross-modal validation, a literature-RAG known-biology
checklist (identity determined by the agent, never told), fault-injection eval, and
novel findings with adversarial Judger review. Full worked examples, both generated on `claude-opus-5`:
[`reports/examples/tenx-cell-ranger-report.md`](reports/examples/tenx-cell-ranger-report.md) (total cost
$14.48) and [`reports/examples/shareseq-multi-cell-lines-report.md`](reports/examples/shareseq-multi-cell-lines-report.md)
(total cost $42.35 — the private dataset's 8-cell-line pool means substantially more
tool-call context per step than the public single-sample dataset; both reports' own
Cost breakdown sections show the real per-step split). The private report's own
fault-injection section caught something genuinely useful beyond the injected faults: a
real RNA-side vs. ATAC-side cell-line label discordance in the dataset itself, found
independently across every one of Opus's four scenario runs including the unmodified
clean control.

This replaced an earlier standalone 4-model (Haiku/Sonnet/Opus/Fable) cost-vs-quality
comparison across the whole eval as the project's primary deliverable — that standalone
table/plot is retired. But the fault-injection section specifically still compares all
4 models head-to-head (the one axis worth comparing even within a single-model report):
the same fault-scenario set is run once per model and rendered as a comparison table
inside each comprehensive report, alongside the report's own model's full per-scenario
answers for qualitative detail. See PROGRESS.md for the full history.

## Datasets used in the example reports

- **Public — 10x Genomics PBMC 10k Multiome.** A real, publicly downloadable same-cell
  RNA+ATAC dataset of healthy peripheral blood mononuclear cells — one unlabeled donor
  sample, no cell-line or treatment-condition metadata at all. `bash
  scripts/download_pbmc_data.sh` fetches it directly from 10x Genomics. This is also the
  dev/demo dataset: small subsamples (500–3000 cells) keep iteration fast.
- **Private — an unpublished multi-cell-line, multi-condition multiome dataset (the
  author's own).** Real same-cell RNA+ATAC data pooling **8 distinct cell lines**, with
  genotype-confirmed ground-truth cell-line labels available for scoring only (e.g.
  cell-line-recovery ARI) — never told to the agent, which has to determine identity from
  the data's own fields and markers like any other run. This dataset has no drug/treatment
  axis (verified by `check_for_condition_groups`, not assumed) — every cell is a control.
  Per this project's data-handling rules, **the data itself never leaves the local
  machine, and its source — file paths, project name, lab/institution name — is never
  named anywhere in this repo.** A real cell-line name the agent itself determines from
  evidence during a run may appear in that run's output (cell-line names aren't treated
  as sensitive, only the source is); see "Stretch phase" below for the full policy.

## Stretch phase: private multi-cell-line data

Beyond the public PBMC dev dataset, this project's full pipeline was also validated end to
end on a private, real multi-cell-line multiome dataset — not a bigger PBMC run, genuinely
different data (pooled cancer cell lines, different QC fields, no ATAC fragments file, no
raw RNA counts layer). Per this project's data-handling rules, **the private data never
leaves the local machine, and its SOURCE (file paths, project/lab/institution names) is
never named or described anywhere in this repo.** A real cell-line name the agent itself
determines from evidence during a run may appear in that run's output — a deliberate,
scoped exception (cell-line names aren't treated as sensitive, only the source is) — but
nothing else identity-bearing is ever committed. See `PROGRESS.md`'s "Stretch" sections
for full detail; every number below is real and independently verified before being
written here.

- **Pipeline validation on real, different data.** RNA-vs-ATAC cross-modal cluster
  agreement: **ARI = 0.725** (vs. PBMC's 0.467 — distinct cell lines separate more cleanly
  than PBMC subtypes, a sensible result). This dataset's direct ground-truth analog to
  "known-biology recovery" — does unsupervised clustering recover the real, genotype-based
  cell-line identity? — gives **ARI = 0.767 (RNA), 0.906 (ATAC)** vs. true label, both far
  above chance.
- **Two new fault types** the single-sample PBMC data couldn't exercise: **cell-line label
  swap** (a demultiplexing-error simulation) degrades cell-line-recovery ARI cleanly and
  monotonically with severity — **0.767 → 0.623 → 0.358 → 0.157** at 0/10/30/50% of cells
  swapped (plotted below). **Mixed samples** (two sequencing libraries merged/mislabeled as
  one) shows a real but modest QC-spread effect — reported honestly as modest, not oversold.
- **4-model fault-detection comparison, run against these two fault types plus shuffled
  RNA-ATAC pairing** (now Section 4 of [`shareseq-multi-cell-lines-report.md`](reports/examples/shareseq-multi-cell-lines-report.md)
  itself, not a standalone comparison file — that standalone artifact was retired in the
  resume-project rebuild in favor of folding the comparison into each comprehensive
  report): all four models detect and correctly diagnose both original faults (2/2), with
  the same convergent
  "false alarm" pattern seen on the public data — every model independently flags a real,
  elevated doublet rate on the unmodified clean data, a genuine QC observation about real
  data, not a fluke. This run also caught and fixed a real bug: the shared system prompt
  hardcoded "this is a PBMC dataset" for every run, including private-data ones — a live eval
  transcript caught a model (Opus) explicitly flagging the resulting "premise mismatch" before
  the fix, and the fix itself measurably changed results (Haiku's detection rate improved once
  it wasn't confused about what dataset it was looking at).
- **chromVAR-style motif deviations, initially wrongly deferred, then corrected**: an early
  scope decision skipped motif analysis for this dataset citing "no peak annotation" — that
  conflated peak *annotation* (nearest-gene mapping, genuinely absent) with peak
  *coordinates* (all that motif scanning actually needs, and already parseable from the same
  `var_names` format as the public data). Corrected once caught: the same `hg38.2bit`/JASPAR
  infrastructure now runs directly against this dataset's ~134K filtered peaks, real result
  879 motifs retained across 5,814 cells, with `tf_motif_correlation` verified working
  end-to-end. Building this also surfaced a real, separate bug in the *shared* public
  pipeline code (an in-place division that silently failed on unsigned-integer counts, masked
  on the public data only because that loader happens to produce float32) — fixed once, now
  more robust for both datasets.
- **Ablation (fixed core only vs. core + agent-chosen analyses)**: the agent's real value-add
  was interpretation and synthesis (combining several QC numbers into one coherent argument
  with concrete follow-ups) and honest tool-choice calibration (proactively searching
  literature, then honestly reporting when the search was inconclusive rather than
  fabricating support) — not new facts beyond the fixed-core tables. A real, useful, honest
  answer to "does agent choice help," not a foregone positive result.
- **Split-half / negative-control novelty checks**: split-half replication **succeeded** —
  two independently-reclustered random halves converged on the same finding as the full
  dataset. The negative control (60% cell-line label swap) **did not cleanly pass**: the
  agent correctly read the fault-corrupted recovery number but confidently explained it as
  real biology rather than a data-integrity artifact — a genuine, disclosed limitation (the
  agent has no way to distinguish "unusual because of biology" from "unusual because the
  ground truth was corrupted" from a downstream summary number alone), not smoothed over.
- **Agent-written analysis code**: a deliberate, bounded exception to this project's "the
  agent never writes new analysis code" MVP rule. Runs locally (never Anthropic's
  server-side code-execution tool, which would require uploading data to Anthropic's
  sandbox — incompatible with "summaries, not raw data, go to the LLM"), restricted to
  pre-aggregated summary statistics only, never raw per-cell data. Building it surfaced a
  real security gap: the initial denylist blocked `open(` but not `pandas`/`numpy`'s own
  file I/O (`pd.read_csv`, `np.load` read arbitrary local files even with `open(` blocked) —
  found by direct testing, fixed by denylisting both libraries' full I/O surface, and
  covered by a dedicated regression test.

## Quickstart

### 1. See the demo (public PBMC data, fastest path)

```bash
conda env create -f environment.yml       # or: pip install -e .[dev]
cp .env.example .env                      # add your ANTHROPIC_API_KEY (never committed)
bash scripts/download_pbmc_data.sh        # ~3.3GB: PBMC 10k Multiome + hg38 2bit genome
make demo                                 # one real agent investigation, ~10s, ~$0.01
```

`download_pbmc_data.sh` fetches the public 10x Genomics PBMC 10k Multiome dataset and
the UCSC hg38 2bit reference genome into `data/raw/` (gitignored — not part of this
repo). It's idempotent, safe to re-run. Note: the **full fixed-core pipeline** (QC,
clustering, gene activity, motif deviations on the raw data) takes **~17-20 minutes on
first run**, dominated by `snapatac2` sorting the whole ~2GB fragments file — a fixed
cost regardless of subsample size. `make demo` sidesteps this by using an
already-committed cache built during development (`data/processed/`, also gitignored —
regenerated automatically, but the *first* real run anywhere still pays the ~17-minute
cost once). This is disclosed, not hidden: see PROGRESS.md throughout for exactly where
this cost shows up and why.

```bash
make test    # full test suite (pytest) — the expensive fixed-core tests also take ~20-30 min
```

Don't want to run anything yourself? The two committed reports in `reports/examples/`
are the finished output of exactly this pipeline — read those directly.

### 2. Generate the full comprehensive report yourself

This is a separate, explicit step (not wired into a Makefile target — it costs real
money across several agent calls and, for the public dataset, one ~17-minute
fault-injection rerun, so it shouldn't run by accident):

```bash
python -c "from multiome_agent.agent.report_generator import generate_report; \
            generate_report('tenx-cell-ranger', model='claude-opus-5')"
```

### 3. Using this on your own data

The fault-injection section is the most expensive part of `generate_report()` and
assumes things arbitrary data may not have (a raw ATAC fragments file for downsampling,
a raw RNA counts layer for doublet injection) — **you don't need it** to get value from
the rest of the pipeline. QC, agent-decided loading, identity discovery, the
literature-RAG known-biology checklist, and adversarially-judged novel findings all work
standalone against any same-cell RNA+ATAC data:

1. Point `config/local_paths.yaml` (gitignored) at your own files, using the same keys
   the private dataset uses — this works for any data split across a separate RNA file
   (all genes), an optional RNA-highly-variable-genes file (used only for
   clustering/embeddings — omit and point it at the same file as the all-genes one if you
   don't have a separate HVG file), and a separate ATAC/peaks file:
   ```yaml
   shareseq_rna_h5ad: /path/to/your_rna_all_genes.h5ad
   shareseq_rna_hvg_h5ad: /path/to/your_rna_highly_variable_genes.h5ad
   shareseq_atac_h5ad: /path/to/your_atac_peaks.h5ad
   ```
   (If your data is instead a single combined 10x-style `.h5` file with both modalities
   together, add a `tenx_matrix_h5` key pointing at your own file instead — this used to
   be hardcoded to one fixed filename, a real limitation now fixed the same way the
   per-modality paths above work:
   ```yaml
   tenx_matrix_h5: /path/to/your_combined_rna_and_atac.h5
   ```
   )
2. Load + let the agent decide the loading strategy from your files' real structure
   (never told which strategy to use):
   ```python
   from multiome_agent.data.dispatch import load_fixed_core_via_agent_decision
   mdata, loader_decision = load_fixed_core_via_agent_decision("shareseq-multi-cell-lines", model="claude-opus-5")
   ```
3. Call just the pieces you want — no fault-injection eval, no 4-model comparison:
   ```python
   from multiome_agent.agent.shareseq_qc_summary import format_shareseq_qc_summary, shareseq_fixed_core_summary
   from multiome_agent.agent.checklist_generator import generate_checklist
   from multiome_agent.agent.novelty import propose_novel_findings, judge_novel_findings
   from multiome_agent.agent.prompts import SHARESEQ_DATASET_CONTEXT

   qc_summary = format_shareseq_qc_summary(shareseq_fixed_core_summary(mdata))
   items, _ = generate_checklist(mdata, qc_summary, SHARESEQ_DATASET_CONTEXT, model="claude-opus-5")
   checklist_summary = "\n".join(f"- ({i.category}) {i.claim} [PMID {i.pmid}]" for i in items)
   findings, _ = propose_novel_findings(mdata, qc_summary, SHARESEQ_DATASET_CONTEXT, checklist_summary, model="claude-opus-5")
   judged = judge_novel_findings(findings, mdata, qc_summary, SHARESEQ_DATASET_CONTEXT, model="claude-opus-5")
   ```
   `SHARESEQ_DATASET_CONTEXT` is generic (no hardcoded specifics about the author's own
   data) — it's exactly the "identity is never told, real cell-line names are OK to
   surface, the data's source never is" policy described above, reusable as-is for
   anyone's own multi-cell-line or multi-condition data. You get back plain Python
   objects (`ChecklistItem`/`JudgedFinding` lists) to inspect, print, or format however
   you like — `generate_report()`'s markdown rendering always bundles all 5 sections
   together, fault-injection included, so it isn't reused here.

## Repo layout

```
src/multiome_agent/
  data/        # loaders (public PBMC + private multi-cell-line), subsampling, caching,
               # dispatch.py (agent-decided loading strategy -> cached fixed-core getters)
  core/        # fixed QC + clustering + gene activity + motif deviations (public);
               # a parallel fixed-core-equivalent pipeline for the private dataset;
               # cross_modal_validation.py, condition_detection.py, identity_metadata.py
  menu/        # agent-selectable analysis-menu tool(s)
  tools/       # PubMed RAG (search + fetch-abstract) / Enrichr wrappers (cached),
               # data_inspection.py (cheap file-structure peek), depmap.py (DepMap
               # Model ID -> cell-line-name resolution), sandboxed code-execution
  agent/       # loop, prompts (identity-discovery principle, near-empty per-dataset
               # context), loader_selector.py, checklist_generator.py, novelty.py
               # (propose + adversarial Judger), report_generator.py (the 5-section report)
  eval/        # fault injection, scoring, grounding checker (public + private-dataset-
               # specific fault types), reused directly by report_generator.py's §4
tests/
reports/examples/   # committed real outputs: the 5-section reports for both datasets --
                     # tenx-cell-ranger-report.md is a full worked example; shareseq-multi-cell-lines-report.md is
                     # aggregate-only per its data-handling rules (see Data-handling above)
scripts/             # data download, demo, plotting
```

The private multi-cell-line dataset's own loader/pipeline/eval modules follow the same
layout conventions but are only usable with a local, gitignored `config/local_paths.yaml`
pointing at that data — see "Stretch phase" above. Nothing about that dataset (paths,
project name, identities) appears anywhere in this repo.

## Limitations

These are the real, specific limitations this project actually found while building and
evaluating it — not generic boilerplate:

- **Single-sample public dev dataset.** The PBMC 10k Multiome dataset has no cell-line or
  condition labels, so its own fault-injection benchmark only has three fault types
  (shuffled pairing, ATAC downsampling, injected doublets). The two remaining CLAUDE.md
  fault types (cell-line/cell-type label swaps, mixed samples) are exercised on the private
  multi-cell-line dataset instead — see "Stretch phase" above — not on PBMC.
- **Small subsample sizes (500-3000 cells)** for fast iteration mean some checklist items
  are noisier than they'd be on the full ~12k-cell dataset — PAX5 motif recovery is a
  genuine, disclosed miss at n=500 (p=0.76), not a scorer bug.
- **Heuristic free-text classifiers**, not exact ground-truth matching, score whether an
  agent "detected" a problem or its diagnosis "matches" the true fault — necessarily
  approximate for open-ended text. Two specific, precisely-diagnosed failure modes are
  documented in PROGRESS.md: a substring false-positive (`"problem"` inside
  `"problematic"`) and the grounding checker flagging legitimate background-knowledge QC
  thresholds (e.g. "FRiP > 0.3 is good") the same way it flags real fabricated numbers.
- **Prompt caching is wired correctly but not currently active** for this account —
  verified directly (a 2400+ token literal, byte-identical system prompt still shows
  zero `cache_creation`/`cache_read` tokens on repeated calls), not assumed.
- **(Resolved in the resume-project rebuild)** ~~No cluster-marker-listing tool~~ —
  `top_cluster_markers`/`top_gene_activity_markers` now expose real per-cluster
  differential-expression results on both modalities, which is what identity discovery
  and the known-biology checklist are actually built on.
- **Injected-doublet scenarios don't recompute ATAC clustering/motifs** for the
  synthetic cells (they have no real underlying fragments) — the eval's QC summary for
  that scenario blends RNA-side numbers honestly but ATAC-side numbers still describe
  only the original population.
- **A genuinely interesting, unresolved finding**: Sonnet, Opus, and Fable all
  independently flag a real statistical property of the clean control (median
  mitochondrial content of 9.9% implies ~half the cells exceed a common ~10% cutoff) —
  correct, sophisticated reasoning that nonetheless counts as a "false alarm" under this
  eval's strict binary ground truth. This exposes a real limitation of "clean_control"
  as a label for real (not synthetically perfect) data, more than it exposes a model
  flaw. See PROGRESS.md Day 2 step 7b.
- **Multi-model comparison ran once per model, no repeats** (Haiku alone has a 3x
  consistency check) — cost-bounded by design; a fuller consistency study across all
  four models is future work.
- **(Resolved in the resume-project rebuild)** ~~Gene activity and motif deviations
  aren't built for the private dataset~~ — motif deviations only ever needed peak
  *coordinates* (parseable from `var_names`) plus the existing hg38.2bit/JASPAR infra,
  not the nearest-gene *annotation* this dataset's peaks genuinely lack (an earlier scope
  decision conflated the two). Gene activity's public-data method (`snap.pp.make_gene_matrix`)
  genuinely can't run here (no fragments file), but the standard ArchR/Signac fallback —
  summing peak counts within each gene's window — can, reusing the SAME cached GENCODE
  gene-coordinate reference the public method uses (gene coordinates are a property of
  genome build, not of which cells were sequenced) rather than a fresh download.
- **The negative-control novelty check did not cleanly pass** on the private dataset: given
  a data-integrity fault (cell-line labels corrupted), the agent read the resulting
  abnormal number correctly but confidently explained it with a biological story rather
  than considering the data itself might be wrong — it has no mechanism to distinguish
  "unusual because of real biology" from "unusual because the ground truth was corrupted"
  from a downstream summary number alone.
- **The agent-written-code sandbox is a restricted namespace plus a string denylist, not a
  formal sandbox.** It stops the concrete gaps found so far (import, dunder/attribute
  escapes, direct and pandas/numpy-mediated file I/O) but isn't a proven boundary against a
  maximally adversarial payload — scoped honestly for a cooperative LLM investigating its
  own pre-aggregated summary data, not for running untrusted third-party code.
- **The private data's three files have no explicit condition/treatment column** — cell
  identity (real, genotype-based) and sequencing-library ID are both present, but library ID
  is the finest-grained grouping available and stands in as the replicate/batch unit for any
  cross-condition-style comparison. A real detection mechanism now exists
  (`check_for_condition_groups`, by column name) and the agent must call it and confirm
  control-only rather than assume it — this dataset genuinely has none, verified, not skipped.
- **ATAC downsampling and doublet injection aren't reused as-is for the private dataset's
  fault-injection section** — a real format constraint (no fragments file; RNA `.X` has no
  raw counts layer), not a design choice. Its fault-injection section uses cell-line-label-swap,
  mixed-samples, and shuffled-RNA-ATAC-pairing instead (4 scenarios total, matching the public
  dataset's scenario count even though the specific mechanisms differ).
- **The private-data ablation and novelty checks ran once, not repeated** — same
  cost-bounded-by-design reasoning as the public multi-model comparison above.

## Status

Day 1 (repo scaffold, data loader, fixed-core pipeline, minimal agent loop), Day 2 (fault
injection, known-biology checklist, grounding checker, report writer, Haiku eval, 4-model
comparison), the Stretch phase (private multi-cell-line data: pipeline validation, two new
fault types, ablation study, split-half/negative-control novelty checks, agent-written
analysis code), and the resume-project rebuild (agent-decided loading strategy, identity
discovery, real-abstract literature RAG, cross-modal cell-type-call validation, a Judger
agent for novel-finding review, a drug/condition-detection branch, and the 5-section
comprehensive report that replaced the 4-model comparison as the primary deliverable) are
all complete. See `PROGRESS.md` for the full build log — every real number in this README
traces to a specific, verified entry there.
