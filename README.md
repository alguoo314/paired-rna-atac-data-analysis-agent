# Multiome QC & Hypothesis Agent

An LLM agent that acts like a computational biologist reviewing a new single-cell RNA + ATAC
dataset: it judges whether the data is good enough to trust, checks whether expected
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
Agent-decided loader
        ↓
Fixed-core analyses
        ↓
Agent loop (planner)
        |
        +----------------------+----------------------+----------------------+----------------------+
        |  check identity &    |  detect conditions   |  search literature    |   run an analysis      |
        |   cross-modal fit    |   (treatment axis,   |   (PubMed: search,    |   (TF-motif corr.,     |
        |  (markers, ARI       |       if any)        | then read abstract)   |   pathway enrich, etc.)|
        +----------------------+----------------------+----------------------+----------------------+
        |
        ↓
Judger
(adversarial review of novel claims)
        ↓
5-section report
           
```


## What it produces

The flagship output is one comprehensive report per dataset, generated end to end by the agent:

1. **Dataset summary** — deterministic RNA/ATAC QC, Leiden clustering with
   RNA-vs-ATAC cluster agreement, gene activity scores, and real chromVAR-style motif deviations.
2. **Gene activity + chromVAR + cross-modal validation** — the textbook ArchR/Signac check: is a
   gene's RNA marker status confirmed by independently-computed ATAC accessibility in its real
   cross-modal partner cluster?
3. **Known-biology checklist**, built in two steps to find previously reported biological claims relevant to the data and confirm against actual data: (1) find claims in the
   literature — `search_pubmed` for candidate titles, then `fetch_pubmed_abstracts` to actually
   read the real abstract before citing anything, since a title match alone never counts; (2)
   verify if the claim holds in this dataset. Aims for up to 3 well-grounded claims tested per category (RNA markers, motifs,
   TF-expression-tracks-motif-accessibility).
4. **Fault-injection eval** — corrupted-data scenarios (shuffled RNA-ATAC pairing, ATAC downsampling, injected doublets, cell-line label swaps, mixed sequencing sample) run across four Claude models (Haiku, Sonnet,
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
  confirmed ground truth available for scoring only — never told to the agent. The data source is not disclosed.

## Quickstart

### See the demo (public PBMC data, fastest path)

```bash
conda env create -f environment.yml       # or: pip install -e .[dev]
cp .env.example .env                      # add your ANTHROPIC_API_KEY (never committed)
make demo                                 # one real agent investigation, ~10-25s, ~$0.02
```

`make demo` runs against a small, already-committed cache (`demo_cache/agent_demo_core.h5mu`,
150 cells, ~46MB — see `agent/demo_fixed_core_cache.py`), so it's genuinely fast on a fresh
clone with no data download needed (a real run: 24s, $0.022). It's illustrative only, at a small
subsample — see `reports/examples/` for the real, full-scale analyses this project's flagship
reports are built from. It asks one fixed default question unless you give it your own:

```bash
make demo QUESTION="does CD3E mark T cells in this dataset?"
```

Want the raw data and the full-scale pipeline (needed for the fresh-report-generation command
below)? That's a separate, larger download, not needed for `make demo` itself:

```bash
bash scripts/download_pbmc_data.sh        # ~3.3GB: PBMC 10k Multiome + hg38 2bit genome
```

The *first* real run of the full fixed-core pipeline on that raw data (QC, clustering, gene
activity, motif deviations, at the real 11,909-cell scale) takes on the order of an hour —
dominated by chromVAR's motif-deviation permutation step at full scale, not the ATAC fragments
sort — a one-time, disclosed cost, never triggered by `make demo` itself.

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

The agent workflow will run through QC, cell identity
discovery, the literature-grounded checklist, and adversarially-judged novel findings.


1. Point `config/local_paths.yaml` (gitignored) at your own files. Separate per-modality files:
   ```yaml
   shareseq_rna_h5ad: /path/to/your_rna_all_genes.h5ad
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
   
      
3. Call just the pieces you want on QC and hypothesis discovery on your data:
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
---

**The pipeline also has a working drug/condition-detection branch — per-arm QC and
condition-specific literature grounding for datasets with a real treatment axis — fully built and
tested. It's not exercised in either example report here: the public PBMC data is control-only,
and the private dataset's real condition arm involves actual drug identities that can't be shown
in a public repo.**



## Roadmap: future analysis-menu additions coming soon

The current analysis menu (gene activity, chromVAR motif deviations, TF-motif self-tracking,
per-cluster differential markers) covers the textbook cross-modal validation but not the more
mechanistic questions a real regulatory-genomics analysis would ask next. Candidates for future
validated wrappers, in roughly the order they'd add the most value:

- **Peak-to-gene links (cis-co-accessibility, à la Cicero/ArchR).** Right now, "distal-enhancer
  regulation" is one of the explanations the agent is *allowed* to reach for when RNA and ATAC
  disagree (see core principle 2), but there's no tool that actually tests it — it can only look
  at gene-body/promoter-proximal accessibility. A real peak-to-gene linkage score would let it
  check whether a DISTAL peak's accessibility (not the gene's own gene-activity score) tracks the
  target gene's expression, closing the gap between "explainable discordance" and "explained."
- **Peak-level differential accessibility between two named groups.** The menu currently only
  offers per-cluster marker DE (RNA) and its gene-activity analog (ATAC) — there's no way to ask
  "which peaks differ between cell line A and cell line B" or "which peaks differ between drug and
  vehicle" directly. This would give the condition-detection branch (already built, not yet
  exercised by either example dataset) a real per-arm result to report, not just per-arm cell
  counts.
- **Regulon inference (TF → target-gene-set, SCENIC+-style).** In practice, novelty proposals in
  this project keep collapsing into "TF X's own RNA tracks its own motif" — a narrow,
  easily-already-known pattern, because that's the richest signal the current toolset can surface.
  A real regulon (does
  a TF's RNA level track the COORDINATED expression of its *predicted target genes*, not just
  itself?) is a structurally different, harder-to-already-know claim, and would give the novelty
  step a genuinely different kind of candidate to propose instead of a rescoped version of the
  same TF-motif-self-correlation check the checklist already runs.
