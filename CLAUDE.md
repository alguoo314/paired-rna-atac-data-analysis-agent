# Project: Multiome QC & Hypothesis Agent

An LLM agent that analyzes paired single-cell RNA + ATAC (multiome) data across cell lines/conditions. It (1) judges data quality, (2) checks which expected RNA signatures, ATAC signatures, and RNA-ATAC relationships are recovered, and (3) proposes candidate novel relationships, grounded in database and literature queries with citations. A fault-injection benchmark measures how well it works.

---

## Core design principles

1. **Deterministic stats, LLM reasoning.** All numbers (QC metrics, differential tests, motif scores, correlations) are computed in Python code. The LLM plans, chooses analyses, queries databases, and interprets. The LLM never generates numbers.
2. **Summaries, not raw data, go to the LLM.** Send top genes, QC tables, enrichment results, etc. Never send raw matrices. (Cost + privacy: the owner will later run this on private data.)
3. **Discordance is expected, not a bug.** RNA and ATAC correlate weakly per gene. The agent should classify RNA-ATAC discordances as *explainable* (poised/primed chromatin, distal-enhancer regulation, housekeeping promoters open everywhere, mRNA stability, missing TF) vs. *surprising*, with evidence. Prefer regulatory-level linkages over naive per-gene correlation:
   - TF expression vs. TF motif accessibility (chromVAR-style) -- `menu/tf_motif_correlation.py`. The narrowest of the three: tests a TF against its OWN motif only.
   - Peak-to-gene links (distal peaks vs. nearby gene expression) -- implemented (`core/peak_to_gene_links.py` + the `peak_to_gene_links` tool): tests whether a DISTAL peak (explicitly outside the gene's own gene-activity window) tracks that SAME gene's own expression, turning "distal-enhancer regulation" from a reachable-but-untested explanation into one the agent can actually check. Deliberately scoped to one gene's own expression on one side -- never a transcription factor's expression standing in for it, a different question.
   - Regulon inference, TF -> predicted target-gene-set (SCENIC+/RcisTarget-lite) -- implemented (`core/regulon_inference.py` + the `regulon_inference` tool): for a TF, finds candidate target genes via real motif + chromatin evidence (the TF's motif present in a target's own promoter, OR in a genuinely distal peak that itself clears `peak_to_gene_links`' own significance bar for that specific peak/gene pair), then reports, PER TARGET GENE, the TF's own RNA correlated against that ONE target's own RNA -- deliberately never an aggregate "regulon activity" score across many genes, because a specific named gene pair (e.g. "SPI1 regulates CD14") is what a literature-checking agent can actually look up and confirm; an averaged module score isn't a citable claim. An earlier, aggregate-score design for this was considered and dropped mid-build for exactly that reason before being redesigned per-target-gene -- see PROGRESS_phase2.md for the full back-and-forth.
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
- ATAC aggregation level: gene activity / motif-level / peak-to-gene links / regulon inference (all implemented) / peak-level differential accessibility between two named groups (**not yet implemented** -- see README roadmap; the condition-detection branch below has no per-arm analytical result yet, only per-arm cell counts)
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
2. **Gene activity, chromVAR motif deviations, and cross-modal validation** -- the textbook ArchR/Signac cross-modal cell-type-call check, run as a systematic sweep over EVERY RNA cluster rather than scraped from whichever genes another section happened to investigate: for each cluster's own top-scoring marker gene, is it a significant RNA marker of that specific cluster AND does its independently-computed ATAC gene-activity confirm elevated accessibility in that cluster's real cross-modal partner (matched by cell overlap, not by clusters coincidentally sharing an integer label)? A cluster whose top-scoring gene doesn't clear the standard significance threshold is skipped, not forced with a non-significant marker -- that's a fine, expected outcome, not an error. This costs $0 in LLM calls (`systematic_cross_modal_sweep` in `core/cross_modal_validation.py` calls the underlying deterministic Python wrappers directly, no agent turn involved) at the cost of real compute time proportional to cluster count.
3. **Known-biology checklist via literature RAG** -- the agent first determines this dataset's real identity from evidence, never told directly (per design principle 6). "From evidence" includes the dataset's own metadata: `check_for_identity_columns` checks `.obs` for a real identity-encoding field (e.g. a genotype-confirmed cell-line column) BEFORE falling back to marker-based inference -- and an indirect encoding counts too: a column of Broad Institute DepMap Model IDs ("ACH-XXXXXX") is resolved to a real cell-line name via `resolve_depmap_id` (a live lookup against EBI's public Cellosaurus database -- DepMap's own portal sits behind an interactive bot-check a plain HTTP client can't pass). No ID->name table is hardcoded anywhere in this project; the agent has to notice the ACH-XXXXXX pattern itself and decide to resolve it. If the dataset pools multiple cell lines, the agent identifies the overall composition then spans the checklist across a handful of different lines chosen at random, rather than concentrating everything on whichever single line is easiest or most already-characterized -- spanning means the SET of claims names different lines (claim 1 about line A, claim 2 about line B), not any single claim hedging generically across several lines at once; each claim stays exactly as specific as if the agent had focused on one line the whole time. The agent then AIMS FOR 3 real, currently-read literature claims per category (RNA marker / motif / peak-to-gene distal link via `peak_to_gene_links` / TF-target-gene regulon link via `regulon_inference` -- the last one specifically a DIFFERENT gene than the TF itself, verified against that tool's per-target `significant` flag, not the TF's own motif) -- a genuine target, not a ceiling. (A 5th category, TF-expression-tracks-its-own-motif-accessibility via `tf_motif_correlation`, was retired from the checklist -- too narrow/easily-already-known a claim next to "regulon_target"'s TF-vs-specific-OTHER-gene claim; the tool itself remains real and usable elsewhere.) settling for fewer requires having actually tried multiple candidates first -- each verified present in this dataset via a real tool call before being recorded. High-impact journals preferred when multiple sources exist for the same claim; substitute a different claim if a candidate gene/motif turns out absent from the data. **The checklist reports real rejections, not just successes**: a specific, literature-backed candidate whose own data-verification tool call clearly contradicts the claim (not significant, wrong sign, or absent from this dataset) is recorded too, rendered in the report as its own "Rejected by data" section separate from "Confirmed by data" -- distinct from a candidate that simply checked out fine but wasn't needed once enough confirmed items were found, which isn't recorded at all. **Once identity is determined from evidence, a real cell-line name may appear in the report** (an explicit, scoped loosening of the anonymization default below) -- but nothing about the data's SOURCE ever may.

**Cost-saving shortlist-first path (`agent/shortlist.py`)**: the expensive part of checklist generation was never the PubMed calls (always free) -- it's the open-ended multi-turn exploration (discovering identity from scratch, brainstorming candidates with no hint, trying up to 4 per category). When a dataset's real identity is ALREADY retrievable for free from its own metadata (`detect_identity_columns`), most of that exploration is redundant. `generate_checklist`'s default behavior now tries, PER CATEGORY: a one-shot, tool-free recall on the cheapest model (Haiku) for well-known candidates -> verify each against THIS dataset's real data via the SAME already-free deterministic tools (zero LLM) -> one targeted literature check per surviving candidate (one search + one cheap single-abstract judgment, not an open-ended brainstorm). Only categories that don't reach 3 confirmed items this way get a SCOPED fallback `run_agent` call covering just those categories (with the already-resolved identity passed in, so it isn't re-derived). If no identity is cheaply available at all (e.g. the public PBMC dataset), this has zero effect -- full original behavior, unchanged. `categories`/`existing_items` parameters let a caller add new categories to an already-generated checklist without re-running or re-paying for categories that haven't changed, which come back completely untouched. Deliberately NOT implemented: provenance tagging of which path (cheap vs. expensive) produced a given item -- the report stays silent on that distinction, by explicit instruction.
4. **Fault-injection eval, compared across all 4 Claude models (Haiku/Sonnet/Opus/Fable)** -- public: shuffled RNA-ATAC pairing, ATAC downsampling, injected doublets, plus a clean control. Private: cell-line label swap, mixed samples, shuffled RNA-ATAC pairing, plus a clean control (ATAC downsampling/doublet injection aren't reusable as-is here -- no fragments file, no raw RNA layer -- a real format constraint, documented in the report's limitations, not silently skipped). The scenario set (including the one real ~17-min snapatac2 downsampling rerun for the public dataset) is built ONCE and reused across all 4 models; the section renders a comparison table (faults detected / correct diagnosis / false alarms / cost per model) plus one model's (the report's own) full per-scenario answers for qualitative detail. Those per-scenario answers don't need a separate "Evidence" paragraph breaking down the supporting tool calls -- the model's direct answer and diagnosis is enough; the comparison table above it is already the grounded, numeric evidence for the section.
5. **Novel findings, adversarial judging, and limitations** -- the proposal step itself tries up to 6 candidates beyond the checklist, running its own literature search (`search_pubmed` + `fetch_pubmed_abstracts`) on EACH one before deciding whether to record it: a candidate the model's own search already shows is well-established is silently self-rejected (no tool call for it), not recorded. Recording stops once 3 candidates survive this self-check, or after 6 distinct candidates are tried, whichever comes first -- ending with fewer than 3 (even zero) is a legitimate, expected outcome, not a failure to force past. **Priority order**: the first 3 attempts must be `regulon_inference`- or `peak_to_gene_links`-based (a TF-target-gene or peak-to-gene claim) -- never `tf_motif_correlation`-based (TF-vs-own-motif) -- only falling back to a TF-vs-own-motif candidate for remaining attempts once those first 3 fail to survive their own literature check. Each recorded finding is then independently reviewed by a separate, explicitly skeptical "Judger" agent call that checks for artifacts (confounds, small samples, non-replication across an independent tool call) AND does its own literature RAG to check prior art -- but FIRST, a cheap prefilter (one free search + one Haiku single-abstract judgment) checks whether the finding is trivially already known; if so it's struck down right there and the expensive Judger call is skipped entirely for it. It's fine -- expected, even -- for every finding to still be struck down at that stage too. A negative control repeats the proposal step against an actually-shuffled (not just a doctored QC-summary string) RNA-ATAC pairing and reports whether it hallucinated a spurious cross-modal finding (0 is well-behaved).

**Drug/condition branch:** both current datasets are control-only. `core/condition_detection.py` detects a drug/treatment/condition axis by `.obs` column name (never by inspecting/guessing at values, and excluding identity-bearing columns) -- the agent must call this and confirm control-only rather than assume it. For a future dataset that does have one, per-arm cell-count QC (`condition_group_qc`) and condition-specific literature RAG (papers either cell-line-agnostic or matching a cell line actually in the dataset) become part of the analysis; this is a real, tested branch (synthetic-column test), just not exercised by either current dataset.

**Cost:** each real run is `checklist_generation`/`novelty_proposal`/`negative_control` (single-model, multi-turn) + `judging` (one call per finding) + `fault_injection` (now 4 models x N scenarios each, the most expensive section by call count though each individual call is cheap) + one small `loader_decision` call. Develop/iterate on Haiku; generate the final committed report on a stronger model, and check the actual total cost before committing to it (see `ReportCost` in the module).

---

## Data

- **Development:** 10x Genomics PBMC ~10k cells Multiome (same-cell RNA+ATAC, well characterized). Subsample to a few thousand cells for fast iteration.
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

### 9. Model comparison -- narrowed to fault-injection only
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



## Working conventions for Claude Code

- Work in small, testable steps; run and verify each before moving on.
- Explain design decisions briefly so the owner can defend them in interviews.
- Ask before downloading large files or making many paid API calls.
- Don't silently change evaluation definitions or ground-truth checklists.
