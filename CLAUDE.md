# Project: Multiome QC & Hypothesis Agent

An LLM agent that analyzes paired single-cell RNA + ATAC (multiome) data across cell lines/conditions. It (1) judges data quality, (2) checks which expected RNA signatures, ATAC signatures, and RNA-ATAC relationships are recovered, and (3) proposes candidate novel relationships, grounded in database and literature queries with citations. A fault-injection benchmark measures how well it works.

Goal: a polished, portfolio-quality GitHub repo built in ~1-2 days (MVP first, stretch goals later). The owner is a comp bio PhD student; she must understand and be able to defend every design decision, so keep code readable and explain non-obvious choices in comments/docs.

**Resume-project rebuild (current primary deliverable):** the flagship output is now one comprehensive 5-section report per dataset (public PBMC, private multi-cell-line) -- see "Comprehensive report format" below. This superseded the earlier plan of a README results table comparing Haiku/Sonnet/Opus/Fable on cost vs. quality across a fixed eval set; that comparison is retired as a deliverable (see "Evaluation" §8).

---

## Core design principles

1. **Deterministic stats, LLM reasoning.** All numbers (QC metrics, differential tests, motif scores, correlations) are computed in Python code. The LLM plans, chooses analyses, queries databases, and interprets. The LLM never generates numbers.
2. **Summaries, not raw data, go to the LLM.** Send top genes, QC tables, enrichment results, etc. Never send raw matrices. (Cost + privacy: the owner will later run this on private data.)
3. **Discordance is expected, not a bug.** RNA and ATAC correlate weakly per gene. The agent should classify RNA-ATAC discordances as *explainable* (poised/primed chromatin, distal-enhancer regulation, housekeeping promoters open everywhere, mRNA stability, missing TF) vs. *surprising*, with evidence. Prefer regulatory-level linkages over naive per-gene correlation:
   - TF expression vs. TF motif accessibility (chromVAR-style)
   - Peak-to-gene links (distal peaks vs. nearby gene expression)
   - Per-modality identity signatures, then check cross-modality agreement
4. **Every claim is grounded.** Each claim in the report cites a computed table and/or a database/literature source, and carries a confidence label. "Novel" claims require a statistical test + a literature check showing it's not already well known.
5. **Everything is logged.** The agent's plan, tool calls, and justifications are saved as a readable reasoning trail.
6. **The agent determines dataset identity from evidence, never told directly.** What tissue/cell type(s), cell line(s), or condition(s) a dataset contains is discovered from real tool results (cluster markers, cross-modal validation), not stated in any system prompt. This applies to loading strategy too: the agent inspects real file structure and decides whether a data source is a combined single-file format or separate per-modality files, rather than a human/script pre-selecting which loader to call. Literature grounding follows the same standard: a title match alone is not evidence -- an abstract must actually be fetched and read before a claim can cite it.

---

## Architecture

```
agent-decided loader -> fixed core analyses -> agent loop (planner) -> report writer
  (inspects real file                              |
   structure; picks               +--------------------+--------------------+
   combined vs. per-modality      |                    |                    |
   loading strategy)      analysis menu          database tools       literature tool
                       (validated wrappers)   (also: condition-axis   (search + fetch-abstract
                                                detection, cross-       RAG, not title-only)
                                                modal validation)
```

The loader decision is itself agent-mediated (see `agent/loader_selector.py` + `tools/data_inspection.py`): a small dedicated agent call inspects cheap file metadata (shape, column names, feature-type values -- never the matrix data) and picks a named loading strategy via a schema-constrained tool call; the actual parsing stays fully deterministic Python underneath. For the real report pipeline this hands off to the CACHED fixed-core getters (`data/dispatch.py`'s `load_fixed_core_via_agent_decision`), not a fresh reload, so gene activity/chromVAR aren't recomputed just to demonstrate the loader decision.

### Fixed core (always runs, gives a stable baseline for QC verdicts and eval)
- RNA QC: genes/UMIs per cell, % mito, doublet scores
- ATAC QC: fragments per cell, TSS enrichment, FRiP, nucleosome signal
- Clustering / cell-line assignment; RNA-vs-ATAC cluster agreement
- Default chromVAR-style motif deviation scores
- Gene activity scores

### Analysis menu (agent chooses, must log justification)
Validated wrappers around established tools with constrained parameters. The agent does NOT write new analysis code in the MVP.
- ATAC aggregation level: gene activity / peak-level differential accessibility / motif-level / peak-to-gene links
- Pseudobulk or metacells vs. per-cell (e.g. pseudobulk when ATAC is sparse)
- Choice of contrasts, subclustering, which database to query next

### Stack
Python-native only (avoid R bridges): scanpy, muon, SnapATAC2, pychromVAR, anthropic SDK, python-dotenv. Verify package availability/APIs before relying on them.

### External knowledge tools (verify each API works; cap output sizes)
Candidates: PubMed E-utilities (literature), Reactome and/or g:Profiler/Enrichr (pathways), STRING (protein interactions), Open Targets (drug targets), JASPAR (motifs), DepMap/CCLE (cell-line expression references), ENCODE cCREs (regulatory elements). Truncate results (e.g. top 10 hits, abstracts only) and cache responses locally. PubMed is two calls, not one: `search_pubmed` (titles, cheap, broad) to find candidates, then `fetch_pubmed_abstracts` (real abstract text) on a short, already-narrowed list -- a title match alone is never sufficient grounding for a citation.

---

## Comprehensive report format (`agent/report_generator.py`)

The flagship deliverable per dataset (public PBMC, private multi-cell-line) is one markdown report with 5 sections, generated by `generate_report(source, model)`:

1. **Dataset & fixed-core summary** -- cell/gene/peak counts, QC medians, cluster counts, RNA-vs-ATAC cluster agreement; for the private dataset, also cell-line count and cell-line-recovery ARI. Also states which loading strategy the agent itself decided from real file structure.
2. **Gene activity, chromVAR motif deviations, and cross-modal validation** -- the textbook ArchR/Signac cross-modal cell-type-call check: for each gene the agent examined elsewhere in the run, is it a significant RNA marker of some cluster AND does its independently-computed ATAC gene-activity confirm elevated accessibility in that cluster's real cross-modal partner (matched by cell overlap, not by clusters coincidentally sharing an integer label)? Built from real `cross_modal_marker_check` calls the OTHER sections' agents already made, not a redundant separate step.
3. **Known-biology checklist via literature RAG** -- the agent first determines this dataset's real identity from evidence, never told directly (per design principle 6). "From evidence" includes the dataset's own metadata: `check_for_identity_columns` checks `.obs` for a real identity-encoding field (e.g. a genotype-confirmed cell-line column) BEFORE falling back to marker-based inference -- and an indirect encoding counts too: a column of Broad Institute DepMap Model IDs ("ACH-XXXXXX") is resolved to a real cell-line name via `resolve_depmap_id` (a live lookup against EBI's public Cellosaurus database -- DepMap's own portal sits behind an interactive bot-check a plain HTTP client can't pass). No ID->name table is hardcoded anywhere in this project; the agent has to notice the ACH-XXXXXX pattern itself and decide to resolve it. If the dataset pools multiple cell lines, the agent identifies the overall composition then focuses the checklist on one clearly-characterized line rather than spreading across all of them. The agent then AIMS FOR 3 real, currently-read literature claims per category (RNA marker / motif / TF-expression-tracks-motif-accessibility) -- a genuine target, not a ceiling; settling for fewer requires having actually tried multiple candidates first -- each verified present in this dataset via a real tool call before being recorded. High-impact journals preferred when multiple sources exist for the same claim; substitute a different claim if a candidate gene/motif turns out absent from the data. **Once identity is determined from evidence, a real cell-line name may appear in the report** (an explicit, scoped loosening of the anonymization default below) -- but nothing about the data's SOURCE ever may.
4. **Fault-injection eval, compared across all 4 Claude models (Haiku/Sonnet/Opus/Fable)** -- public: shuffled RNA-ATAC pairing, ATAC downsampling, injected doublets, plus a clean control. Private: cell-line label swap, mixed samples, shuffled RNA-ATAC pairing, plus a clean control (ATAC downsampling/doublet injection aren't reusable as-is here -- no fragments file, no raw RNA layer -- a real format constraint, documented in the report's limitations, not silently skipped). The scenario set (including the one real ~17-min snapatac2 downsampling rerun for the public dataset) is built ONCE and reused across all 4 models; the section renders a comparison table (faults detected / correct diagnosis / false alarms / cost per model) plus one model's (the report's own) full per-scenario answers for qualitative detail.
5. **Novel findings, adversarial judging, and limitations** -- up to 3 candidate novel findings beyond the checklist, each independently reviewed by a separate, explicitly skeptical "Judger" agent call that checks for artifacts (confounds, small samples, non-replication across an independent tool call) AND does its own literature RAG to check prior art. It's fine -- expected, even -- for every finding to be struck down. A negative control repeats the proposal step against an actually-shuffled (not just a doctored QC-summary string) RNA-ATAC pairing and reports whether it hallucinated a spurious cross-modal finding (0 is well-behaved).

**Drug/condition branch:** both current datasets are control-only. `core/condition_detection.py` detects a drug/treatment/condition axis by `.obs` column name (never by inspecting/guessing at values, and excluding identity-bearing columns) -- the agent must call this and confirm control-only rather than assume it. For a future dataset that does have one, per-arm cell-count QC (`condition_group_qc`) and condition-specific literature RAG (papers either cell-line-agnostic or matching a cell line actually in the dataset) become part of the analysis; this is a real, tested branch (synthetic-column test), just not exercised by either current dataset.

**Cost:** each real run is `checklist_generation`/`novelty_proposal`/`negative_control` (single-model, multi-turn) + `judging` (one call per finding) + `fault_injection` (now 4 models x N scenarios each, the most expensive section by call count though each individual call is cheap) + one small `loader_decision` call. Develop/iterate on Haiku; generate the final committed report on a stronger model, and check the actual total cost before committing to it (see `ReportCost` in the module).

---

## Data

- **Development:** 10x Genomics PBMC ~10k cells Multiome (same-cell RNA+ATAC, well characterized). Subsample to a few thousand cells for fast iteration.
- **Showcase (multi-cell-line):** 8 prostate cell lines multiome (RWPE1, RWPE2, PrEC, BPH1, DU145, PC3, 22Rv1, LNCaP; bioRxiv 2024, "Identify Regulatory eQTLs by Multiome Sequencing in Prostate Single Cells"). **Verify processed data is actually downloadable on GEO before committing.**
- **Multi-condition option:** CAT-ATAC (10x Multiome + CRISPR guides, cancer cell lines under drug/genetic perturbation, dasatinib), GEO GSE288996.
- **Private multi-cell-line, multi-condition data (owner's own):** already paired to the same cells.
  Three files (ATAC peaks, RNA all-genes, RNA highly-variable-genes-only), paths configured
  locally via `config/local_paths.yaml`'s `shareseq_atac_h5ad`/`shareseq_rna_h5ad`/
  `shareseq_rna_hvg_h5ad` keys (gitignored -- see `src/multiome_agent/config.py`). Never write
  the real absolute paths into this file or any other committed file (see "Data handling
  rules" below -- this instruction file is part of the repo too).
  - Use the **all-genes** RNA file for markers, differential expression, TF expression, and known-biology checks (many canonical markers and TFs may not be HVGs). The HVG file is fine for clustering/embeddings.
  - First step: inspect both files (obs/var columns, cell-line and condition labels, layers raw vs. normalized, peak coordinates format, whether cell barcodes match exactly and in the same order) and summarize before any analysis.
  - **Data handling rules:**
    - Never copy these files or any derived cell-level data (per-cell matrices, embeddings, barcodes) into the repo.
    - Aggregate outputs and findings (reports, summary tables, plots, eval results) CAN be committed and used as the agent showcase.
    - Never name or describe the data SOURCE anywhere in the repo (no file paths, project names, lab/institution names, or file names). Refer to the source generically, e.g. "an unpublished multi-cell-line, multi-condition multiome dataset".
    - **Cell-line identity is a scoped exception (resume-project rebuild decision):** once the agent determines a real cell-line name from evidence in its own run (not by being told), that name may appear in agent output/reports -- cell-line names themselves aren't treated as sensitive, only the source (file paths, project/lab/institution names) is. Older fault-injection ground-truth code that anonymizes cell-line/library identities to letters (A, B, C...) predates this decision and still anonymizes by construction; that's a stricter, still-safe default for ground-truth data specifically, not a contradiction to resolve.
    - Load paths from a gitignored config (e.g. `.env` or `config/local_paths.yaml`), not hardcoded, and before committing check that no path or source name leaks into code, reports, logs, or plot titles.
    - Only summary statistics go to the LLM API, as with all data.
- Prefer same-cell multiome. Datasets where scRNA and scATAC were profiled separately are less suitable for per-cell RNA-ATAC links.

---

## Evaluation (key selling point of the repo)

### 1. Analysis code validation
Unit tests on synthetic data with planted signals; compare core outputs on PBMC against an established workflow (markers, motif scores broadly agree).

### 2. Fault-injection benchmark
Corrupt copies of clean datasets:
- Cell-line / cell-type label swaps
- Shuffled RNA-ATAC barcode pairing (looks fine per modality, breaks joint structure)
- Heavy downsampling of one modality
- Mixed samples
- Injected doublets
- Optional severity levels (e.g. 5% / 20% / 50% labels swapped)

**Always include clean, uncorrupted controls** to measure false alarms.
Score separately: detected a problem? correct cause diagnosed? clean data left alone?

### 3. Known-biology recovery
**Current mechanism (resume-project rebuild):** the checklist is generated fresh each report run via literature RAG against whatever identity the agent itself determines from evidence (see "Comprehensive report format" §3) -- not scored against one fixed pre-written ground-truth file. This is a deliberate change from the original MVP plan below, made because a static PBMC-only / prostate-only checklist doesn't generalize to "whatever dataset gets pointed at this pipeline," which the identity-discovery design principle now requires anyway.
Original MVP plan (still a valid, complementary evaluation methodology -- write expected answers BEFORE looking at agent output, store as a checklist file, score recall + count incorrect "expected signature" claims -- just not what the comprehensive report itself does):
- PBMC: RNA markers CD14/LYZ (monocytes), MS4A1 (B), CD3E (T), NKG7 (NK); motifs SPI1/PU.1 and CEBP (myeloid), PAX5/EBF1 (B), TBX21/EOMES (NK); TF expression should track its motif accessibility.
- Prostate: AR, KLK3, NKX3-1, FOXA1 activity in LNCaP/22Rv1; absent in AR-negative PC3/DU145.

### 4. Grounding checks (automated)
- Every number in the report matches computed tables.
- Every cited PMID / database ID exists. Owner spot-checks that sources actually support the claims.

### 5. Novel findings
**Current mechanism:** up to 3 candidate findings, each adversarially reviewed by a separate "Judger" agent call (artifact check via independent tool calls + its own literature RAG for prior art) -- see "Comprehensive report format" §5. Negative control: propose-novel-findings run again against an actually-shuffled RNA-ATAC pairing; 0 findings (or none that survive scrutiny) is well-behaved, a hallucinated cross-modal finding from noise is not.
Not yet built: split-half replication across random cell halves, cross-dataset replication, owner's expert rating (plausible / trivial / wrong) as a separate, independent check beyond the Judger's own verdict.

### 6. Consistency
Repeat key evals 3-5 times; report spread, not best run.

### 7. Ablation
Fixed core only vs. core + agent-chosen analyses. Does agent choice actually help?

### 8. Model comparison -- narrowed to fault-injection only (resume-project rebuild decision)
The original plan was a standalone README results table comparing Haiku/Sonnet/Opus/Fable across the WHOLE eval (faults, known-biology recall, citations, refusals), plus a cost-vs-accuracy scatter plot, as its own deliverable. That standalone table/plot is retired -- the comprehensive report (see "Comprehensive report format" above) is the primary deliverable, generated once on one model. But the fault-injection SECTION of that report specifically still compares all 4 models (per explicit direction: this is the one axis worth comparing head-to-head even within the single-model report), reusing the SAME scenario set across all 4 -- so `eval/fault_injection.py`, `eval/scoring.py`, `eval/grounding.py`, and their private-data equivalents remain real, tested, and load-bearing, not vestigial. Removed: the multi-model comparison-table-building code covering the OTHER eval axes and its output artifacts (`model_comparison.md`, `eval_results_*.json`).

---

## Cost & safety conventions

- Model is a single config setting: `AGENT_MODEL` env var, default `claude-haiku-4-5`. Develop on Haiku; demo/eval on stronger models. Do a few full runs on the demo model before final eval (model-specific bugs).
- API key lives ONLY in project `.env` (`ANTHROPIC_API_KEY=...`), loaded with python-dotenv. `.env` is in `.gitignore`. Never export it globally in the shell (Claude Code would then bill the API key instead of the owner's Max subscription).
- Max-turns limit on the agent loop.
- Prompt caching for the long fixed system prompt.
- Cache LLM and database responses locally so reruns of downstream steps are cheap.
- Log token usage and cost per run (needed for the results table anyway).

---

## Repo deliverables (what makes it look good)

- README: one-paragraph pitch, architecture diagram, "results at a glance" pointing to the two comprehensive reports, quickstart, limitations section.
- Two polished comprehensive reports committed in the repo, one per dataset (readers see output without running anything).
- One-command setup (`pip install -e .` or conda env file) and one-command demo (`make demo` or single script).
- Tests for analysis code and tool wrappers.
- Clear module layout, e.g.:
```
src/multiome_agent/
  data/        # loaders, subsampling, dispatch.py (agent-decided loader)
  core/        # fixed QC + baseline analyses, cross_modal_validation.py, condition_detection.py
  menu/        # agent-selectable analysis wrappers
  tools/       # database + literature wrappers (with caching), data_inspection.py
  agent/       # loop, prompts, loader_selector.py, checklist_generator.py,
               # novelty.py (propose + Judger), report_generator.py (5-section report)
  eval/        # fault injection, scoring, grounding (single-model, reused by report_generator)
tests/
reports/examples/
```

---

## Suggested build order

**Day 1**
1. Repo scaffold, env, config, `.env` handling, logging.
2. Data loader + subsampled PBMC multiome.
3. Fixed core analyses + unit tests; sanity-check against known PBMC biology.
4. Minimal agent loop with 2-3 tools (one analysis-menu item, one database, literature), max-turns guard, reasoning log. Test on Haiku.

**Day 2**
5. Fault injector + clean controls + scoring harness.
6. Known-biology checklist + grounding checkers + report writer.
7. Run eval on Haiku, then model comparison on reduced set.
8. README, example report, plots, limitations.

**Stretch (later)**
Multi-cell-line showcase dataset with my private data, ablation, split-half/negative-control novelty checks, agent-written analysis code.

**Resume-project rebuild (after Stretch)**
Replaced the model-comparison deliverable with the 5-section comprehensive report described above: agent-driven loader selection, identity discovery (never told directly), literature RAG with real fetched abstracts (not title matching), cross-modal cell-type-call validation via gene activity, a Judger agent for novel-finding review, and a drug/condition-detection branch. See "Comprehensive report format" for what actually ships.

---

## Working conventions for Claude Code

- Work in small, testable steps; run and verify each before moving on.
- Explain design decisions briefly so the owner can defend them in interviews.
- Ask before downloading large files or making many paid API calls.
- Don't silently change evaluation definitions or ground-truth checklists.
