# Multiome QC & Hypothesis Agent

An LLM agent that reviews a new single-cell **RNA + ATAC (multiome)** dataset the way a computational biologist would. It:

1. **judges whether the data can be trusted** (QC and cross-modal consistency),
2. **checks whether known biology actually shows up** (literature-grounded checklist), and
3. **proposes new candidate findings** for a human to review.

It then writes everything up as one complete report per dataset, end to end, with no human running analysis steps in between.

## Design principles

- **The LLM never produces a number.** Every QC metric, cluster assignment, motif score and correlation comes from validated Python tools. The LLM only decides which analyses to run and interprets the results.
- **No hints about the sample.** The agent is never told the cell type, cell line or condition. It has to infer identity from evidence, the way a scientist reading an unlabeled dataset would.
- **Citations must be read, not matched.** A literature claim counts only after the agent fetches and reads the real PubMed abstract. A title match is never enough.
- **Every claim is traceable.** A grounding checker verifies that each cited number and PMID traces back to a real tool call from the same run.
- **Tested against failure.** A fault-injection benchmark checks that the agent notices planted data problems, and a negative control on fully shuffled data checks that it doesn't hallucinate relationships.
- **Findings face an adversary.** Each proposed novel finding is re-investigated by an independent Judger agent looking for artifacts and prior art.

## Architecture


```
Agent-decided loader ........ inspects file structure, picks a loading strategy
        │
        ▼
Fixed-core analyses ......... QC · clustering · gene activity · chromVAR
        │
        ▼
┌─► Agent loop (planner) ──── picks one tool per step:
│       ├─ Check identity & cross-modal fit   (markers, ARI)
│       ├─ Detect conditions                  (treatment axis, if any)
│       ├─ Search literature                  (PubMed search → read abstract)
│       └─ Run an analysis                    (regulons, peak-gene links, enrichment …)
│               │
└───────────────┘ results return to the planner until it's done
        │
        ▼
Judger ...................... adversarial review of novel claims
        │
        ▼
5-section report
```

## Output: a 5-section report

| # | Section | What it contains |
|---|---|---|
| 1 | **Dataset summary** | Deterministic RNA/ATAC QC, Leiden clustering with RNA-vs-ATAC cluster agreement, gene activity scores, and chromVAR motif deviations |
| 2 | **Cross-modal validation** | The standard ArchR/Signac check: is a gene's RNA marker status confirmed by independently computed ATAC accessibility in its matched cross-modal cluster? |
| 3 | **Known-biology checklist** | Claims from the literature, tested against this dataset (see below) |
| 4 | **Fault-injection eval** | Corrupted-data scenarios run across four Claude models (Haiku, Sonnet, Opus, Fable) side by side |
| 5 | **Novel findings** | Candidate discoveries beyond the checklist, each cross-examined by the Judger, plus a negative control on fully shuffled data |

**How the checklist is built.** The agent first finds candidate claims with `search_pubmed`, then reads the actual abstract with `fetch_pubmed_abstracts` before citing anything. It then tests whether each claim holds in this dataset, aiming for up to 3 well-grounded claims per category:

- RNA markers
- Motifs
- Distal peak-to-gene links
- TF → target-gene regulon links

*Design note:* "a TF's expression tracks its own motif" was retired as a checklist category in favor of regulon links. It is a narrow and usually already-known pattern, whereas a TF tracking a specific *other* gene it is predicted to regulate is a more informative test. The TF-vs-own-motif analysis is still available as a tool.

**Fault-injection scenarios:** shuffled RNA-ATAC pairing, ATAC downsampling, injected doublets, cell-line label swaps, and a mixed sequencing sample.

📄 Two full worked reports are in [`reports/examples/`](reports/examples/).

## Analysis methods

| Analysis | Method / package | Question it answers |
|---|---|---|
| RNA QC | scanpy (`pp.calculate_qc_metrics`) + Scrublet | Genes, UMIs and % mito per cell; predicted doublet score |
| ATAC QC | SnapATAC2 (fragments import) | Fragments per cell, FRiP, TSS enrichment, nucleosome signal |
| RNA clustering | scanpy (HVGs → PCA → Leiden) | Cell groupings from the transcriptome alone |
| ATAC clustering | muon (TF-IDF → LSI, depth component dropped) → Leiden | Cell groupings from chromatin accessibility alone |
| Cross-modal cluster agreement | Adjusted Rand Index (RNA vs. ATAC Leiden labels) | Do the two modalities independently agree on cell groupings? |
| Gene activity | SnapATAC2 `make_gene_matrix` (fragments-based); peak-summation fallback when no fragments file exists | ArchR/Signac-style proxy: accessibility near a gene's promoter and body |
| chromVAR motif deviations | pychromvar + JASPAR 2024 CORE (879 vertebrate motifs) + MOODS scanning | Per cell and motif: is chromatin at the motif's binding sites more or less open than a GC/accessibility-matched background? (TF-activity proxy) |
| Cross-modal marker validation | Custom cell-overlap-matched RNA/ATAC cluster pairing | Is a gene's RNA marker status confirmed by elevated ATAC gene activity in the matched cluster? |
| TF-motif self-tracking | Spearman (TF RNA vs. its own chromVAR deviation) | Does a TF's transcript level track its own motif activity? |
| Peak-to-gene links | Spearman (distal peak accessibility vs. gene RNA), BH-corrected and effect-size-gated | Does one specific distal peak, outside the gene's own activity window, track that gene's expression? (Cicero/ArchR-style) |
| Regulon inference (TF → target) | JASPAR motif in a target's promoter or significant distal peak, then Spearman (TF RNA vs. target RNA) | Does a TF's RNA track a specific gene it is predicted to regulate? (SCENIC+/RcisTarget-lite) |
| Pathway enrichment | Enrichr (GO Biological Process) | What is a gene set, such as cluster markers, collectively involved in? |
| Literature grounding | PubMed E-utilities (`search_pubmed` + `fetch_pubmed_abstracts`) | Is a claim reported in a real abstract that was actually read? |
| Cell-line identity resolution | EBI Cellosaurus via DepMap Model ID lookup | Resolves a coded `ACH-XXXXXX` ID in the metadata to a real cell-line name |
| Fault-injection eval | Synthetic corruption scenarios across 4 Claude models | Does the agent notice planted data problems, and how does that vary by model? |
| Adversarial judging | Independent Judger agent | Is a proposed finding likely an artifact, or already known? |

## Datasets in the example reports

- **Public: 10x Genomics PBMC 10k Multiome.** Publicly available same-cell RNA + ATAC data from healthy peripheral blood mononuclear cells. One unlabeled donor with no cell-line or condition metadata, so the agent has to reconstruct cell-type identity from markers alone.
- **Private: an unpublished multi-cell-line multiome dataset (the author's own).** Same-cell RNA + ATAC data pooling several cell lines. Genotype-confirmed ground truth is used only for scoring and is never shown to the agent. The data source is not disclosed.

## Quickstart

### 1. Set up

```bash
conda env create -f environment.yml       # or: pip install -e .[dev]
cp .env.example .env                      # add your ANTHROPIC_API_KEY (never committed)
```

### 2. Run the demo (no data download)

```bash
make demo
```

The demo runs on a small cache committed to the repo, so it works on a fresh clone. A typical run takes about **14 seconds** and costs about **$0.009**.

By default it asks one fixed question, chosen to showcase the `regulon_inference` tool:

> Does GATA3 regulate any specific target gene in this dataset, not just its own motif?

Ask your own with `QUESTION=`:

```bash
make demo QUESTION="does CD3E mark T cells in this dataset?"
make demo QUESTION="does SPI1 regulate CD14 in this dataset?"
make demo QUESTION="is there a distal enhancer for CD3E?"   # exercises peak_to_gene_links
```

- **Data:** `demo_cache/agent_demo_core.h5mu` (150 cells, ~85 MB), built by `agent/demo_fixed_core_cache.py`
- **Tools:** the same 13-tool toolset as the full pipeline

> **Note:** The demo is illustrative only, since it runs on a small subsample. The full-scale analyses are in [`reports/examples/`](reports/examples/). If you'd rather not run anything, those reports are the finished output of this exact pipeline.

### 3. Run tests

```bash
make test
```

### 4. Full-scale pipeline (optional)

Download the raw data (not needed for the demo):

```bash
bash scripts/download_pbmc_data.sh        # ~3.3 GB: PBMC 10k Multiome + hg38 2bit genome
```

The first full run of the fixed-core pipeline (QC, clustering, gene activity, TF-gene linkage on all 11,909 cells) takes roughly an hour, mostly chromVAR's motif-deviation permutation step. This is a one-time cost and is never triggered by `make demo`.

Generate a full report (this makes several paid agent calls):

```bash
python -c "from multiome_agent.agent.report_generator import generate_report; \
           generate_report('tenx-cell-ranger', model='claude-opus-5')"
```

## Using your own data

Two input layouts are supported, and you don't need to say which one you have. The agent inspects your files' structure and picks the loading strategy itself:

- a single combined 10x-style `.h5` (RNA and ATAC together), or
- separate per-modality files (an RNA `.h5ad` and an ATAC `.h5ad`).

**1. Point `config/local_paths.yaml` (gitignored) at your files.**

Separate per-modality files:

```yaml
shareseq_rna_h5ad: /path/to/your_rna_all_genes.h5ad
shareseq_atac_h5ad: /path/to/your_atac_peaks.h5ad
```

Or one combined file:

```yaml
tenx_matrix_h5: /path/to/your_combined_rna_and_atac.h5
```

**2. Load the data.**

```python
from multiome_agent.data.dispatch import load_fixed_core_from_local_config

mdata, loader_decision = load_fixed_core_from_local_config(model="claude-opus-5")
```

**3. Run QC, the checklist, and novel-finding discovery.**

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

## Treatment / condition branch

The pipeline also includes a fully built and tested drug/condition-detection branch, with per-arm QC and condition-specific literature grounding for datasets that have a real treatment axis. It isn't shown in the example reports: the public PBMC data is control-only, and the private dataset's treatment arm involves drug identities that can't be published.

## Build log
The development history is split into two logs, both written as the work happened (every bug, wrong turn and fix, with evidence):

- [`PROGRESS.md`](PROGRESS.md) covers the core pipeline, built before peak-to-gene links and regulon inference were added.
- [`PROGRESS_phase2.md`](PROGRESS_phase2.md) covers those two analyses and the changes that came with them.