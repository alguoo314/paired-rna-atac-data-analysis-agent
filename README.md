# Multiome QC & Hypothesis Agent

An LLM agent that acts like a computational biologist reviewing a new single-cell RNA + ATAC
(multiome) dataset: it judges whether the data is good enough to trust, checks whether expected
biology actually shows up in it, and proposes new candidate findings for a human to review — then
writes all of that up as one complete report per dataset, end to end, with no human running
analysis steps in between.

Every number in that report — QC metrics, cluster assignments, motif scores, correlations — comes
from validated Python code; the LLM only plans which analyses to run and interprets what comes
back, it never invents a number. It's also never told in advance what cell type, cell line, or
condition the data contains — it has to work that out from evidence, the way a scientist reading
an unlabeled dataset would. A fault-injection benchmark and a negative control check that it
actually notices planted data problems and doesn't hallucinate a relationship out of thin air.

[`PROGRESS.md`](PROGRESS.md) is the full build log, written as the work happened — every bug
found, every wrong turn, every fix, with receipts. Worth a read if you like seeing how the sausage
actually gets made.

## Architecture

```
agent-decided loader -> fixed-core analyses -> agent loop (planner) -> 5-section report
                                                         |
           +----------------------+----------------------+----------------------+----------------------+
           |   run an analysis    |  search literature   |   check identity &   |  detect conditions   |
           |   (TF-motif corr.,   |   (PubMed: search,   |   cross-modal fit    |   (treatment axis,   |
           |   pathway enrich.)   | then read abstract)  | (markers, ARI check) |       if any)        |
           +----------------------+----------------------+----------------------+----------------------+
```

Loading strategy and dataset identity are both agent-decided, never hardcoded or told: a dedicated
tool call inspects real file structure to pick how to load the data, and a separate one works out
what the data actually *is* from cluster markers and metadata — no human tells it either. The agent
loop itself is a hand-rolled tool-use loop (eleven tools across the four categories above, a
max-turns guard, a full reasoning trail per run), not the SDK's built-in tool runner.

## What it produces

The flagship output is one comprehensive report per dataset, generated end to end by the agent:

1. **Dataset & fixed-core summary** — deterministic RNA/ATAC QC, Leiden clustering with
   RNA-vs-ATAC cluster agreement, gene activity scores, and real chromVAR-style motif deviations
   (full JASPAR vertebrate motif set, actual hg38 sequence), plus how the agent decided to load
   and identify the data.
2. **Gene activity + chromVAR + cross-modal validation** — the textbook ArchR/Signac check: is a
   gene's RNA marker status confirmed by independently-computed ATAC accessibility in its real
   cross-modal partner cluster?
3. **Known-biology checklist**, built in two steps for every candidate claim: (1) find it in the
   literature — `search_pubmed` for candidate titles, then `fetch_pubmed_abstracts` to actually
   read the real abstract before citing anything, since a title match alone never counts; (2)
   verify it in this dataset — a real tool call (`tf_motif_correlation`, `cross_modal_marker_check`,
   `top_cluster_markers`) confirms the claim actually holds here, not just in the paper. A claim
   only makes the checklist if both steps pass; if step 2 fails, the agent tries a different
   candidate instead of forcing it. Aims for 3 verified claims per category (RNA markers, motifs,
   TF-expression-tracks-motif-accessibility).
4. **Fault-injection eval** — corrupted-data scenarios (shuffled RNA-ATAC pairing, downsampling,
   injected doublets, dataset-specific label faults) run across four Claude models (Haiku, Sonnet,
   Opus, Fable) side by side, with a grounding checker verifying every cited number/PMID traces
   back to a real tool call from that same run.
5. **Novel findings** — candidate discoveries beyond the checklist, each cross-examined by an
   adversarial Judger agent for artifacts and prior art, plus a negative control against
   fully-shuffled data to check for hallucination.

Two full worked examples are committed in [`reports/examples/`](reports/examples/) — real output,
not mockups.

## Datasets used in the example reports

- **Public — 10x Genomics PBMC 10k Multiome.** A real, publicly downloadable same-cell RNA+ATAC
  dataset of healthy peripheral blood mononuclear cells. One unlabeled donor sample, no cell-line
  or condition metadata — the agent has to reconstruct cell-type identity from markers alone.
- **Private — an unpublished multi-cell-line, multi-condition multiome dataset (the author's
  own).** Real same-cell RNA+ATAC data pooling several distinct cell lines, with genotype-
  confirmed ground truth available for scoring only — never told to the agent. Per this project's
  data-handling rules, **the data itself never leaves the local machine, and its source (file
  paths, project name, lab/institution name) is never named anywhere in this repo.** A real
  cell-line name the agent itself determines from evidence during a run may appear in that run's
  output — cell-line names aren't treated as sensitive, only the source is.

## Quickstart

### See the demo (public PBMC data, fastest path)

```bash
conda env create -f environment.yml       # or: pip install -e .[dev]
cp .env.example .env                      # add your ANTHROPIC_API_KEY (never committed)
bash scripts/download_pbmc_data.sh        # ~3.3GB: PBMC 10k Multiome + hg38 2bit genome
make demo                                 # one real agent investigation, ~10s, ~$0.01
```

`make demo` runs against an already-committed cache, so it's fast. It asks one fixed default
question unless you give it your own:

```bash
make demo QUESTION="does CD3E mark T cells in this dataset?"
```

The *first* real run of the full fixed-core pipeline anywhere (QC, clustering, gene activity,
motif deviations) takes longer — dominated by sorting the raw ATAC fragments file — a one-time,
disclosed cost, not hidden.

```bash
make test    # full test suite
```

Don't want to run anything yourself? The two committed reports in `reports/examples/` are the
finished output of exactly this pipeline — just read those.

Want the full comprehensive report, generated fresh? That's one line, run explicitly since it
costs real money across several agent calls:

```bash
python -c "from multiome_agent.agent.report_generator import generate_report; \
            generate_report('tenx-cell-ranger', model='claude-opus-5')"
```

### Using this on your own data

Both a single combined file (one 10x-style `.h5` with RNA and ATAC together) and separate
per-modality files (an RNA `.h5ad` and an ATAC `.h5ad`) are supported — you don't tell the agent
which one you have. It inspects your files' real structure and decides for itself which loading
strategy applies (see "Architecture" above).

The fault-injection section is the most expensive part of `generate_report()`, and it assumes
things arbitrary data may not have (a raw ATAC fragments file, a raw RNA counts layer). You don't
need it to get value from the rest of the pipeline — QC, agent-decided loading, identity
discovery, the literature-grounded checklist, and adversarially-judged novel findings all work
standalone against any same-cell RNA+ATAC data:

1. Point `config/local_paths.yaml` (gitignored) at your own files. Separate per-modality files:
   ```yaml
   shareseq_rna_h5ad: /path/to/your_rna_all_genes.h5ad
   shareseq_rna_hvg_h5ad: /path/to/your_rna_highly_variable_genes.h5ad   # optional
   shareseq_atac_h5ad: /path/to/your_atac_peaks.h5ad
   ```
   Or a single combined 10x-style `.h5` file with both modalities together:
   ```yaml
   tenx_matrix_h5: /path/to/your_combined_rna_and_atac.h5
   ```
2. Load, letting the agent decide the loading strategy from your files' real structure:
   ```python
   from multiome_agent.data.dispatch import load_fixed_core_from_local_config
   mdata, loader_decision = load_fixed_core_from_local_config(model="claude-opus-5")
   ```
   No dataset-name argument here either — this auto-detects which config slot from step 1 you
   actually filled in and loads accordingly.
3. Call just the pieces you want — no fault-injection eval, no 4-model comparison:
   ```python
   from multiome_agent.agent.shareseq_qc_summary import format_shareseq_qc_summary, shareseq_fixed_core_summary
   from multiome_agent.agent.checklist_generator import generate_checklist
   from multiome_agent.agent.novelty import propose_novel_findings, judge_novel_findings

   qc_summary = format_shareseq_qc_summary(shareseq_fixed_core_summary(mdata))
   items, _ = generate_checklist(mdata, qc_summary, model="claude-opus-5")
   checklist_summary = "\n".join(f"- ({i.category}) {i.claim} [PMID {i.pmid}]" for i in items)
   findings, _ = propose_novel_findings(mdata, qc_summary, checklist_summary=checklist_summary, model="claude-opus-5")
   judged = judge_novel_findings(findings, mdata, qc_summary, model="claude-opus-5")
   ```
   No dataset-context argument to pass or import: with none given, the agent is told not to
   assume an identity, and once it works one out itself, it may say a real name out loud (like a
   cell line) but never where the data came from — the right default for anyone's own
   multi-cell-line or multi-condition dataset. You get back plain Python objects
   (`ChecklistItem`/`JudgedFinding` lists) to inspect or format however you like.
4. Fault-injection eval is the most expensive piece (each scenario is its own full agent run) —
   most people won't need it standalone, but if you do want to test it against your own data:
   ```python
   from multiome_agent.eval.shareseq_eval_harness import build_shareseq_scenarios, run_shareseq_model_eval

   scenarios = build_shareseq_scenarios(mdata)  # clean control + cell-line-swap + mixed-samples + shuffled-pairing
   result = run_shareseq_model_eval("claude-opus-5", scenarios, mdata)
   print(f"${result['total_cost_usd']:.4f} total")
   for r in result["scenario_results"]:
       print(r["fault_type"], "detected=", r["detected"], "diagnosis_matches=", r["diagnosis_matches"])
   ```
   `build_shareseq_scenarios` assumes your data has the same `cell_line_name`/`library` columns
   the cell-line-swap and mixed-samples faults need — if it doesn't, `shuffle_rna_atac_pairing`
   (also used to build the shuffled-pairing scenario) still works standalone against any
   same-cell RNA+ATAC data.

---

The pipeline also has a working drug/condition-detection branch — per-arm QC and
condition-specific literature grounding for datasets with a real treatment axis — fully built and
tested. It's not exercised in either example report here: the public PBMC data is control-only,
and the private dataset's real condition arm involves actual drug identities that can't be shown
in a public repo.
