# Multiome QC & Hypothesis Agent report (private)

*Model: `claude-opus-5` · total cost: $42.3547*

## 1. Dataset & fixed-core summary

**Data source:** private multi-cell-line multiome data. Loading strategy was determined by the agent itself from real file structure, not pre-specified: it decided `separate_per_modality_files` -- "No file has a feature_types column mixing "Gene Expression" and "Peaks"; instead the modalities are split across files sharing the same 5,814 obs: two gene-expression matrices (19,129 and 3,697 vars, with gene-specific var columns like mt/ribo/hb/highly_variable/binomial_deviance) and one chromatin-accessibility matrix (1,578,279 vars with only an index column, plus ATAC-specific obs fields such as frip, n_fragment, tsse, frac_mito)."

5814 cells pooled from 8 distinct cell lines, 19129 genes, 133743 ATAC peaks (after feature filtering). RNA: 14 Leiden clusters, median 4258 genes/cell, 9296 UMIs/cell, 7.0% mito, 10.6% doublets (scrublet-based call). ATAC: 10 Leiden clusters, median 8484 fragments/cell, median FRiP 0.62, median TSS enrichment 14.6. RNA-ATAC cross-modal cluster agreement (ARI): 0.725. Recovery of true (genotype-confirmed) cell-line identity from unsupervised clustering (ARI vs. ground truth): RNA=0.767, ATAC=0.906. chromVAR-style motif accessibility deviations are available for this dataset -- use the tf_motif_correlation tool to check a specific transcription factor's expression against its own motif's accessibility.

---

## 2. Gene activity, chromVAR motif deviations, and cross-modal validation

Gene activity scores and chromVAR-style motif deviations are part of the fixed-core pipeline (computed once, cached, reused here). The table below is the textbook ArchR/Signac cross-modal cell-type-call validation: for each gene checked while investigating this dataset (identity discovery, known-biology checklist, novel-finding review), is it a significant RNA marker of some cluster, AND does its independently-computed ATAC gene-activity confirm elevated accessibility in that cluster's real cross-modal partner?

| Gene | RNA marker of cluster | Matched ATAC cluster | Gene-activity confirms |
|---|---|---|---|
| ESR1 | 4 | 4 | True |
| GREB1 | 4 | 4 | True |
| GRHL2 | 4 | 4 | True |
| PGR | 4 | 4 | True |
| PRLR | 4 | 4 | True |
| TRPS1 | 4 | 4 | True |

6/6 checked markers cross-validate between modalities. Per CLAUDE.md's design principle, discordance here is expected, not a bug -- gene activity is a noisier, indirect accessibility proxy, and RNA-ATAC agreement is not assumed to be perfect.

---

## 3. Known-biology checklist (literature RAG, identity determined by the agent)

**Identity determination:** **This is not a tissue sample — it's a pooled panel of 8 human cancer cell lines**, and the dataset's own `.obs` says so directly. Both the RNA and ATAC modalities carry a `cell_line_name` field plus DepMap Model IDs (`ACH-XXXXXX`), with a companion QC flag (`may_have_wrong_cell_line_label_*` = "no" for all 5,814 cells):

| Cell line | DepMap ID | Cells (ATAC label) | Lineage (from `atac_lineage`) |
|---|---|---|---|
| NCI-H838 | ACH-000416 | 1,626 | Lung |
| HCC-44 | ACH-000667 | 817 | Lung |
| SJSA1 | ACH-000748 | 773 | Bone |
| **T-47D** | **ACH-000147** | **697** | **Breast** |
| OVTOKO | ACH-000663 | 627 | Ovary/Fallopian tube |
| LN-229 | ACH-000595 | 579 | CNS/Brain |
| YKG1 | ACH-000570 | 531 | CNS/Brain |
| A-673 | ACH-000052 | 164 | Bone |...

Each item below required the agent to: search PubMed, fetch and actually read a real abstract (not just a title), and verify with a real tool call that the gene/motif is usable in this dataset before recording it.

### RNA markers

- **ESR1**: ESR1 (ERα) is expressed in T-47D, an ERα-positive luminal breast cancer line; here ESR1 is a top RNA marker of the T-47D cluster (RNA cluster 4, logFC 5.01, padj 8.3e-245) with ATAC gene-activity confirmation. (PMID 27888136, *The Journal of steroid biochemistry and molecular biology*, 2017; data-presence verified: True)
- **PGR**: T-47D is progesterone-receptor (PGR)-positive; here PGR is a top RNA marker of the T-47D cluster (RNA cluster 4, logFC 5.32, padj 6.3e-242) and is cross-modally confirmed by ATAC gene activity. (PMID 35617163, *PloS one*, 2022; data-presence verified: True)
- **GREB1**: GREB1, a direct ER-regulated early-response gene, correlates with ER-positive phenotype across breast cancer cell lines and so should be enriched in T-47D; here GREB1 is a top RNA marker of the T-47D cluster (RNA cluster 4, logFC 4.85, padj 2.5e-219) with ATAC gene-activity confirmation. (PMID 11103799, *Cancer research*, 2000; data-presence verified: True)

### Motifs

- **MA0148.5.FOXA1**: The FOXA1 forkhead motif should show elevated chromatin accessibility in T-47D, since FOXA1 is a pioneer factor that exclusively initiates chromatin opening at its own genomic binding sites and is essential for growth of breast cancers; the FOXA1 motif (MA0148.5) is present in this dataset's chromVAR deviations. (PMID 41808995, *bioRxiv : the preprint server for biology*, 2026; data-presence verified: True)
- **MA0037.5.Gata3**: The GATA3 motif should show elevated accessibility in luminal ER+ cells such as T-47D: GATA3 is a luminal-defining transcription factor whose motif accessibility marks the luminal chromatin landscape and is lost upon luminal-to-basal reprogramming; the GATA3 motif (MA0037.5) is present in this dataset's chromVAR deviations. (PMID 35021081, *Cell reports*, 2022; data-presence verified: True)
- **MA0112.4.ESR1**: Estrogen response elements (the ESR1 motif) should be in accessible chromatin in T-47D: in T-47D cells specifically, chromatin accessibility at EREs is predictive of ER genomic binding and regulatory activity (with FOXA1 and GATA3 as the T-47D-specific predictive TFs); the ESR1 motif (MA0112.4) is present in this dataset's chromVAR deviations. (PMID 42469015, *Genome research*, 2026; data-presence verified: True)

### TF-expression-tracks-motif-accessibility

- **FOXA1**: FOXA1 abundance should track accessibility at its own motif in breast cancer cells: acute pharmacological degradation of FOXA1 shows it exclusively initiates chromatin opening at its own genomic binding sites, so FOXA1 levels causally set FOXA1-motif accessibility. In this dataset FOXA1 RNA vs FOXA1 motif (MA0148.5) chromVAR deviation gives Spearman rho = 0.108, p = 2.0e-16. (PMID 41808995, *bioRxiv : the preprint server for biology*, 2026; data-presence verified: True)
- **NR2F2**: NR2F2 expression level should track accessibility at NR2F2-bound/motif-containing ERalpha regulatory elements in ER+ breast cancer: perturbation of NR2F2 expression decreases ERalpha DNA binding and chromatin opening, and NR2F2 co-binds with FOXA1/GATA3 at 85% of ERalpha sites. In this dataset NR2F2 RNA vs its motif (MA1111.2) gives Spearman rho = -0.018, p = 0.18 (not significant). (PMID 31588232, *Theranostics*, 2019; data-presence verified: True)
- **GATA3**: GATA3, a known pioneer factor in ER+ breast cancer that co-binds FOXA1/NR2F2 and renders ERE-containing sites pre-accessible, should have its expression track its own motif's accessibility. In this dataset GATA3 RNA vs GATA3 motif (MA0037.5) gives Spearman rho = -0.079, p = 2.0e-09 -- i.e. weakly negative rather than positive, so the expected tracking is not observed. (PMID 31588232, *Theranostics*, 2019; data-presence verified: True)


---

## 4. Fault-injection eval (compared across models)

| Model | Faults detected | Correct diagnosis | False alarms | Cost/run |
|---|---|---|---|---|
| claude-haiku-4-5 | 1/2 | 2/2 | 0/1 | $0.0358 |
| claude-sonnet-5 | 1/2 | 2/2 | 0/1 | $0.6495 |
| claude-opus-5 | 2/2 | 2/2 | 0/1 | $6.5808 |
| claude-fable-5 | 2/2 | 2/2 | 0/1 | $2.7294 |

With the classifier fixed, **no model's clean-control answer mentions a cell-line-identity or
doublet problem anymore.** Opus's one remaining false alarm is a genuinely new, different
observation (a suspicious RNA genes-per-UMI ratio, discussed below) -- not the cell-type/doublet
false-positive class the last two fixes targeted, and a real, defensible thing to flag on its own
merits.

The fixed classifier also surfaced something more interesting than a bug: **"correct diagnosis"
(2/2 for every model) can now legitimately exceed "faults detected" (1/2 for Haiku and Sonnet)**,
because the two checks measure different things -- "correct diagnosis" only checks whether the
answer's vocabulary touches the right topic (e.g. "cell line," "identity," "recovery" for the
label-swap scenario), regardless of what the model concluded about it. Reading the actual text:
Haiku's label-swap answer explicitly says the ARI collapse "is not a failure of the data or
clustering -- it reflects the underlying biology," and Sonnet's shuffled-pairing answer calls the
cross-modal ARI collapse "a bit counterintuitive" but offers "a plausible, non-alarming
explanation." Both models are discussing the right numbers and the right topic, and both
explicitly explain the injected fault away as benign rather than flagging it -- a real, disclosed
miss on 2 of the 8 fault-scenario attempts, not a scoring artifact.

### Per-scenario answers (claude-opus-5)

**clean_control** (Unmodified private multi-cell-line data.)

## Answer

**The data is mostly high quality, but one number needs caution: the RNA library-complexity ratio.** The ATAC side is clean; the RNA side has a likely counts-matrix problem plus a milder over-clustering issue.

**What the dataset is** (determined from its own metadata, not assumed): 5,814 cells pooled from 8 genotype-confirmed cell lines — NCI-H838 (1,626), HCC-44 (817), SJSA1 (773), T-47D (697), OVTOKO (627), LN-229 (579), YKG1 (531), A-673 (164) — spread over 3 batches, with **no drug/treatment/condition axis** (control-only).

### The specific problem
`get_qc_summary` reports **median 4,258 genes/cell on median 9,296 UMIs/cell**. That implies ~2.2 counts per detected gene, i.e. essentially every one of >4,000 genes detected as a near-singleton. Real droplet libraries at ~9,300 UMIs typically detect ~2,500–3,500 genes, because a large fraction of molecules is consumed by a few very high-expressing genes (mito/ribosomal), leaving fewer molecules to spread across the tail. A ratio of 0.46 genes-per-UMI is not physically impossible but is far flatter than any real count distribution should be.

Most likely causes, in order: (1) the matrix is **not raw counts** — it looks like a normalized/denoised/imputed layer (imputation inflates detected-gene counts while total "UMIs" stays at library depth); (2) genes and UMIs were **computed on different matrices** (e.g. gene count pre-filtering, UMI sum post-filtering of the 19,129-gene set); (3) counts were downsampled after detection was tallied. This matters because any depth-sensitive downstream step (HVG selection, DE, doublet scoring) would be operating on a mis-specified count model.

### Secondary observations
- **RNA recovers ground-truth identity worse than ATAC** (ARI 0.767 vs 0.906). This is backwards from the usual expectation — RNA normally separates cell lines more cleanly than sparser ATAC. Part of this is benign: there are 14 RNA clusters for 8 true lines, and ARI penalizes splitting. I confirmed splitting rather than mixing: RNA clusters 12 and 13 share top markers (TGFBR3, DENND5B, and MUC16 in cluster 12, consistent with the ovarian line OVTOKO), i.e. one line fragmented across two clusters. So the RNA ARI gap is probably mostly a Leiden-resolution choice — but it's also consistent with the count-matrix issue above, and I can't fully separate the two from the QC summary alone.
- **Median 7.0% mito** is on the high side for cultured lines but within normal acceptance limits.
- **A-673 has only 164 cells** — enough to call identity, thin for any per-line quantitative claim.
- The QC tool's own description mentions per-cell doublet scores, but **no doublet rate was returned** in the summary, so I can't assess doublet contamination.
- ATAC QC is solidly good: median 8,484 fragments/cell, FRiP 0.62, TSS enrichment 14.6, 133,743 peaks.
- The metadata carries `may_have_wrong_cell_line_label_*` flags for both modalities, and **all 5,814 cells are flagged "no"** — no suspected mislabels.
- Biology is recoverable despite the above: `cross_modal_marker_check` on MDM2 shows it is a significant RNA marker of cluster 1 with independent ATAC gene-activity confirmation in the matched ATAC cluster 1, and cluster 1's top markers (MDM2 lfc 3.92, COL1A1, COL11A1, COL1A2) match the MDM2-amplified osteosarcoma line SJSA1.

### Recommendation
Before trusting any RNA quantitative result, verify that the RNA layer used for QC is the raw integer count matrix and that genes/cell and UMIs/cell were computed on the same matrix after the same filtering. ATAC-side analyses can proceed as-is.

### Evidence

| Claim | Source | Confidence |
|---|---|---|
| 8 genotype-confirmed lines, names and per-line cell counts; no mislabel flags | `check_for_identity_columns` (`cell_line_name`, `may_have_wrong_cell_line_label_*` all "no") | **High** |
| Control-only, no treatment axis | `check_for_condition_groups` (`is_control_only: true`, no columns found) | **High** |
| 4,258 genes vs 9,296 UMIs median; 7.0% mito; ATAC 8,484 frags / FRiP 0.62 / TSS 14.6; ARIs 0.725 / RNA 0.767 / ATAC 0.906; 14 RNA vs 10 ATAC clusters | `get_qc_summary` | **High** (numbers), **Medium-high** (that the genes:UMI ratio is anomalous — grounded in the tool's numbers plus general expectations for droplet libraries, which is **background knowledge**, not a tool output) |
| Which specific cause (imputed matrix vs. mismatched filtering) | inference only | **Low** — can't discriminate without the raw matrix |
| RNA cluster over-splitting (12/13 same line) | `top_cluster_markers` clusters 12 & 13 (shared TGFBR3, DENND5B; MUC16 lfc 2.83 in 12) | **Medium** |
| Biology still recoverable; cluster 1 = SJSA1 | `top_cluster_markers` cluster 1 (MDM2 lfc 3.92, p_adj 0.0) + `cross_modal_marker_check` MDM2 (RNA marker + ATAC gene-activity confirms) | **High** |
| Doublet rate unassessed | not returned by `get_qc_summary` | **High** (stated as a gap, not a finding) |

Note: two `search_pubmed` queries on scRNA-seq gene-detection/UMI benchmarking returned **zero hits**, so no citation backs the complexity expectation — it is labeled background knowledge above.

**cell_line_label_swap** (1744/5814 cells' cell-line label reassigned to a different real cell line.)

## Answer

**Something does look wrong — but it's a clustering/labeling-recovery problem, not a raw-data-quality problem.**

**First, what this sample is (found in the data, not assumed):** the dataset's own `.obs` carries a genotype-confirmed identity field. It is a pool of **8 cancer cell lines**, 5,814 cells: NCI-H838 (1,626), HCC-44 (817), SJSA1 (773), T-47D (697), OVTOKO (627), LN-229 (579), YKG1 (531), A-673 (164). The coded DepMap IDs resolve consistently (e.g. ACH-000416 → NCI-H838, lung adenocarcinoma). There is **no drug/condition axis** — `check_for_condition_groups` returned none, so this is control-only.

**The raw QC is good.** Median 4,258 genes and 9,296 UMIs per cell; ATAC median 8,484 fragments, FRiP 0.62, TSS enrichment 14.6 — all comfortably above standard multiome pass thresholds. RNA–ATAC cluster agreement is reasonable (ARI 0.725).

**The specific problem:** unsupervised clustering recovers the true cell-line identity very poorly — **ARI 0.363 (RNA) and 0.414 (ATAC) against the genotype-confirmed labels**. Pooled, genetically distinct cell lines are the *easiest* possible clustering problem (they differ by copy number, driver mutations and whole expression programs), and should separate at ARI ≈ 0.9+. Getting ~0.4 means the clustering is not tracking cell line.

Two contributing causes, one of which I could demonstrate directly:

1. **Over-clustering / fragmentation (demonstrated).** There are 14 RNA Leiden clusters for 8 lines. RNA clusters 12 and 13 share their top markers (TGFBR3, DENND5B, MCTP2, plus MUC16 in cluster 12) — they are two sub-clusters of one underlying line, not two distinct lines. Splitting single lines across multiple clusters mechanically depresses ARI.
2. **Probable residual mixing as well.** Fragmentation alone doesn't usually push ARI as low as 0.36, so some clusters likely also pool cells from more than one line — plausibly a batch effect (the metadata carries a 3-level batch field, 3,186 / 1,324 / 1,304 cells) or a shared stress/quality axis. I could not cross-tabulate cluster × line with the available tools, so this part is inference, not a measured result.

Notably, the clustering *can* resolve lines when it works: RNA cluster 1 is unambiguously SJSA1 — top marker **MDM2 (LFC 3.92, p_adj ≈ 0)** alongside COL1A1/COL11A1/COL1A2, matching MDM2-amplified osteosarcoma (background knowledge). So the signal is present in the data; the clustering step is failing to exploit it.

**Two smaller caveats:** median mitochondrial fraction is 7.0%, on the high side for cultured lines (tool result + background knowledge); and the A-673 arm has only 164 cells, too few to support confident line-level claims. Also, the QC tool advertises doublet metrics but **returned no doublet number**, so I can't assess doublet burden — I'm not going to guess one.

**Recommendation:** re-run clustering at lower resolution and/or with batch integration before using cluster labels for anything, and benchmark against the genotype labels already in the metadata (which are the trustworthy identity source here).

### Evidence

| Claim | Source | Confidence |
|---|---|---|
| 8 genotype-confirmed cell lines, counts as listed; control-only | `check_for_identity_columns` (`cell_line_name`, `cell_lines`), `check_for_condition_groups` (`is_control_only: true`) | **High** |
| ACH-000416 = NCI-H838, lung adenocarcinoma | `resolve_depmap_id` (Cellosaurus CVCL_1594) | **High** |
| Raw RNA/ATAC QC passes (4,258 genes, 9,296 UMIs, FRiP 0.62, TSS 14.6) | `get_qc_summary` | **High** |
| Identity recovery is poor: ARI 0.363 (RNA) / 0.414 (ATAC) vs ground truth | `get_qc_summary` | **High** |
| Over-clustering: 14 RNA clusters for 8 lines; clusters 12 & 13 share TGFBR3/DENND5B/MCTP2 markers → same line split in two | `list_clusters`, `top_cluster_markers` (12, 13) | **High** |
| Cluster 1 = SJSA1 (MDM2 LFC 3.92, p_adj ≈ 0, + collagen program) | `top_cluster_markers` (1) + background knowledge on MDM2-amplified osteosarcoma | **Medium-high** |
| Residual batch/mixing as an additional ARI driver | Inference from batch field + ARI gap; **not** directly measured | **Low-medium** |
| 7.0% median mito is somewhat elevated for cell lines | `get_qc_summary` + background knowledge | **Medium** |
| Doublet burden | **Not returned by the tool** — not assessed | — |

**shuffled_rna_atac_pairing** (2907/5814 cells have their ATAC profile swapped with another cell's.)

## Answer

**Something does look wrong — but it is not cell-quality, it's the cross-modal (RNA↔ATAC) integration.**

**What the dataset is** (discovered from its own `.obs`, not assumed): 5,814 cells pooled from 8 genotype-labelled cancer cell lines — NCI-H838 (1,626), HCC-44 (817), SJSA1 (773), T-47D (697), OVTOKO (627), LN-229 (579), YKG1 (531), A-673 (164). No drug/treatment axis exists (`check_for_condition_groups` → control-only). Every cell is flagged `may_have_wrong_cell_line_label... = "no"` in both modalities, so the label QC is clean.

**The per-modality data are good.** ATAC: median 8,484 fragments/cell, FRiP 0.62, TSS enrichment 14.6 — all comfortably above standard thresholds. RNA: median 4,258 genes and 9,296 UMIs/cell. Biology within each modality is coherent: RNA cluster 4 is unambiguously T-47D (ESR1 LFC 5.01, PGR 5.32, GREB1 4.85, TRPS1 4.15, PRLR 6.46, all adj-p < 1e-218).

**The specific problem: the QC numbers are internally inconsistent about cross-modal agreement.**
- RNA clustering recovers true cell-line identity at ARI 0.767; ATAC clustering recovers it at ARI 0.906.
- Yet RNA-vs-ATAC cluster agreement is ARI 0.190.

Two clusterings that each track the same 8-way ground truth that closely cannot disagree with *each other* that badly. The 14-vs-10 cluster granularity mismatch does cost some ARI — but the ground-truth ARIs already absorb that penalty (RNA's 14 clusters only cost it ~0.23 against truth), so granularity does not plausibly explain a drop to 0.19. This points at the per-cell RNA↔ATAC barcode correspondence / integration step, not at the cells.

**Corroborating (weaker) signal:** TF expression vs. its own motif's chromVAR deviation is essentially null, and twice the wrong sign, for TFs that should be sharply lineage-restricted in this pool — GATA3 rho = −0.079, FLI1 rho = −0.051, FOXA1 rho = +0.108. With T-47D contributing ~12% of cells, a real per-cell pairing should give a clearly positive GATA3 correlation. I weight this only moderately: motif-based TF-activity inference is genuinely noisy and method-dependent (PMID 39441876 benchmarks exactly this and finds results vary substantially by method), and small |rho| between modalities is normal in general.

**Counter-evidence I have to report honestly:** the pairing is *not* fully scrambled. `cross_modal_marker_check` matched partners by cell overlap and confirmed ESR1 and GATA3 (RNA cluster 4 ↔ ATAC cluster 4) and FLI1 (RNA 3 ↔ ATAC 8), 3/3 with independent gene-activity confirmation. A completely shuffled barcode mapping would have broken those. So this reads as *partially degraded or mis-computed* cross-modal correspondence, not total corruption.

**Two smaller flags:**
- Median 7.0% mitochondrial reads. For a nuclei-based multiome prep this is high (background knowledge: snRNA-seq is typically ≲5%, often ~1%) and suggests ambient/cytoplasmic carryover. Not disqualifying, but worth checking.
- A-673 has only 164 cells (2.8%). Too few for confident line-specific claims; restrict detailed work to the well-powered lines (NCI-H838, HCC-44, SJSA1, T-47D).
- The QC summary did not return a doublet rate despite offering one — I can't comment on doublets.

## Evidence

| Claim | Source | Confidence |
|---|---|---|
| 8 named cell lines, counts as listed; no condition axis; no mislabelled cells | `check_for_identity_columns`, `check_for_condition_groups` | **High** — read directly from real metadata fields |
| Per-cell QC (fragments 8,484; FRiP 0.62; TSS 14.6; 4,258 genes; 9,296 UMIs) is acceptable | `get_qc_summary` | **High** |
| Cross-modal ARI 0.190 is inconsistent with ground-truth ARIs 0.767 / 0.906 → integration/pairing defect | `get_qc_summary` (all three numbers) + reasoning about ARI behaviour (background knowledge) | **Medium** — the numbers are solid; the *inference* that pairing is at fault is the interpretive step |
| Near-null / wrong-sign TF-motif correlations (GATA3 −0.079, FLI1 −0.051, FOXA1 +0.108) | `tf_motif_correlation` ×3 | **High** for the values; **Low–Medium** as evidence of a defect, since motif-activity inference is intrinsically noisy (PMID 39441876) |
| Pairing is not *fully* broken | `cross_modal_marker_check` on ESR1, GATA3, FLI1 — 3/3 confirmed | **Medium-High** |
| Cluster 4 = T-47D | `top_cluster_markers("4")` + metadata | **High** |
| 7.0% mito is high for a nuclei prep | `get_qc_summary` value + background knowledge on snRNA-seq mito content | **Medium** |

**Recommended next step:** re-derive the RNA and ATAC cell barcode index alignment before trusting any joint RNA–ATAC analysis (peak-gene links, motif-expression correlation, multi-omic embedding). Per-modality clustering and cell-line-level marker analysis can be used as-is.

## 5. Novel findings, adversarial judging, and limitations
(Some findings omitted due to private dataset)

### Finding (confidence: high)

GRHL2 — not FOXA1, GATA3 or ESR1 — is the transcription factor whose own RNA best predicts its own motif's chromVAR accessibility in this dataset, by a ~3.5-fold margin, and GRHL2 is simultaneously a strong T-47D-cluster marker confirmed in both modalities. This nominates GRHL2 (a PR co-regulator at distal enhancers) rather than the canonical luminal pioneers as the dominant accessibility-tracking factor in this pooled panel.

**Evidence:** tf_motif_correlation: GRHL2 RNA vs MA1105.3.GRHL2 deviation Spearman rho = 0.373, p = 4.6e-191 — versus FOXA1 0.108 (p=2.0e-16, established), ESR1 0.082 (p=3.9e-10), GATA3 -0.079 (p=2.0e-09, established), TFAP2C 0.054 (p=4.4e-05), SPDEF 0.038 (p=3.9e-03), PGR 0.030 (p=0.021), XBP1 -0.016 (p=0.21, n.s.), ELF3 0.003 (p=0.84, n.s.), KLF5 -0.040 (p=2.5e-03). Other epithelial-restricted TFs (ELF3, KLF5, SPDEF) are near zero, so the effect is not a generic 'epithelial TF has high dynamic range in a pooled panel' artifact. top_cluster_markers(4): GRHL2 logFC 3.35, padj 2.2e-203; cross_modal_marker_check(GRHL2): RNA marker of cluster 4, matched ATAC cluster 4, gene-activity confirms. Dataset context: check_for_identity_columns shows 8 genotype-labeled lines, T-47D = ACH-000147, 697 cells; check_for_condition_groups: control-only. Literature: PMID 41843622 (abstract read) reports GRHL2 interacts with PR and is co-recruited with PR to distal enhancers in hormone-responsive breast cancer cells — consistent with, but not equal to, this accessibility-tracking claim.

**Judger verdict: struck_down**

- Likely artifact: True -- The rho is reproducible (GRHL2 0.373, p=4.6e-191) but is a cross-cell-line identity statistic, not a within-cell-type regulatory readout. check_for_identity_columns shows 8 genotype-labeled lines across 5 lineages with T-47D contributing only 697/5814 cells, so ~88% of cells driving the GRHL2 correlation are non-breast; check_for_condition_groups confirms control-only. Two unexcluded boring explanations: (1) prevalence-capped dynamic range - ESR1/PGR/GATA3/FOXA1 are essentially confined to the single 697-cell line (top_cluster_markers(4): ESR1 logFC 5.01, PGR 5.32, GATA3 4.73) whereas GRHL2's lower logFC 3.35 indicates expression shared beyond that line; a Spearman correlation against a near-binary predictor is structurally bounded by the 'on'-group prevalence (background statistical knowledge), giving a mechanical multi-fold advantage with no biological difference. (2) Motif-family degeneracy, not GRHL2 specialness - the cited controls fail because ELF3 (reproduced rho 0.003, p=0.84) and KLF5 (-0.040) have degenerate ETS/GC-box motifs shared with ubiquitously expressed family members, and GATA3's NEGATIVE rho (-0.079) shows GATA motif deviation is driven by other GATA factors in non-breast lines; the GRHL motif uniquely lacks a ubiquitously expressed family surrogate here (GRHL1 rho 0.062; GRHL3 has no motif in this dataset). Critically, the claimed ~3.5-fold margin collapses once other lines' lineage TFs are tested, which the original evidence omitted: PAX8 rho = 0.214 (p=4.4e-61, ovarian OVTOKO, family-unique motif) and RUNX2 = 0.129, so the true margin is ~1.7x and the generic pattern 'lineage-restricted TF + family-unique motif => high rho' reproduces in unrelated lineages (also TEAD1 0.088, POU3F2 0.096, TFAP2C 0.054, SOX9 -0.044, FLI1 -0.051, TP63 -0.006). The cross-modal leg is true but non-discriminating: cross_modal_marker_check confirms GRHL2 (RNA cluster 4 / ATAC cluster 4), but ESR1 and CDH1 confirm identically and top_gene_activity_markers(4) lists FOXA1 (1.78) and GATA3 (1.89) alongside GRHL2 (1.93), so the canonical pioneers pass the same test equally well. Data quality is not the issue (FRiP 0.62, TSS enrichment 14.6, ATAC-vs-truth ARI 0.906); the inferential design is.
- Already known: True -- Abstracts actually read, not just titles. PMID 29867222 (Nat Genet 2018) causally establishes that Grainy head binding sites determine epithelial enhancer accessibility, that Grh loss/ectopic expression causes loss/gain of DNA accessibility, that human GRHL1/GRHL2/GRHL3 function similarly, and concludes Grh binding is 'necessary and sufficient for the opening of epithelial enhancers' - a tight GRHL2-expression-to-GRHL2-motif-accessibility coupling is the direct predicted consequence of this established mechanism. PMID 32974388 (review) states explicitly that GRHL factors act as pioneer factors establishing a cell-type-specific accessible chromatin landscape exclusive to epithelial transcription. The 'rather than the canonical luminal pioneers' framing is also not new: PMID 31644911 already designates GRHL2 a lineage-determining factor collaborating with FOXA1 in ER+/endocrine-resistant breast cancer, and PMID 36036613 shows GRHL2 is pre-bound at chromatin and required for maximal ER recruitment at enhancers. Thus GRHL2 being the factor whose expression best tracks its own motif accessibility is an expected restatement of well-established pioneer-factor biology, not a novel nomination.

### Finding (confidence: medium)

Both GATA-family factors assayed here show NEGATIVE self-motif tracking, and the stronger of the two is TRPS1 — the atypical GATA repressor, which is the single highest-ranked RNA marker of the T-47D cluster (ranked above ESR1, PGR and GATA3) and is cross-modally confirmed. This offers a concrete, mechanistic candidate explanation for the otherwise anomalous negative GATA3 RNA-vs-GATA3-motif correlation: GATA-element accessibility in this luminal line is co-occupied by a NuRD-recruiting repressor rather than being a pure GATA3-activation readout.

**Evidence:** tf_motif_correlation: TRPS1 RNA vs MA1970.2.TRPS1 deviation Spearman rho = -0.148, p = 7.4e-30 (more strongly negative than the already-established GATA3 rho = -0.079, p = 2.0e-09). top_cluster_markers(4): TRPS1 is the #1-ranked marker of the T-47D RNA cluster, logFC 4.15, padj 6.9e-275, above ESR1 (5.01, 8.3e-245), PGR (5.32, 6.3e-242) and GATA3 (4.73, 4.8e-177) by DE score. cross_modal_marker_check(TRPS1): RNA marker of cluster 4, matched ATAC cluster 4, gene activity confirms. Literature (abstracts read): PMID 30563971 shows TRPS1 is an atypical GATA factor that recognizes GATA elements and represses transcription by recruiting CHD4/NuRD(MTA2), including enhancer decommissioning at TP63; PMID 19759027 confirms TRPS1 binds a GATA consensus site directly to repress a target promoter. Caveat: the correlations are computed across all 5814 cells of an 8-line pool, so cross-line variance contributes; the tools available cannot restrict the correlation to T-47D cells alone.

**Judger verdict: struck_down**

- Likely artifact: True -- The stated premise is factually wrong in this dataset, and the inference is confounded by the pooled design.

(1) PREMISE FALSE. The claim says "both GATA-family factors assayed here show NEGATIVE self-motif tracking." I found five assayable GATA-family factors, and the majority are POSITIVE: TRPS1 -0.148 (p=7.4e-30), GATA3 -0.079 (p=2.0e-09), but GATA2 +0.134 (p=1.2e-24), GATA6 +0.139 (p=2.5e-26), GATA4 +0.045 (p=6.4e-04) (tf_motif_correlation). GATA1/GATA5 are absent from the RNA var_names. There is no GATA-family-wide negative self-motif tracking.

(2) MECHANISM FALSIFIED BY (1). GATA2/3/4/6 and TRPS1 all read essentially the same WGATAA consensus (background knowledge), so their chromVAR deviation scores are largely redundant measurements of the same accessibility feature. If GATA-element accessibility were globally confounded by a NuRD-recruiting co-occupying repressor, every GATA-motif deviation would behave the same way. Instead the sign flips according to which TF's RNA it is — which is the signature of a between-cell-line contrast, not of chromatin co-occupancy.

(3) POOLED-LINE CONFOUND IS THE WHOLE EFFECT. check_for_identity_columns shows an 8-line pool (NCI-H838 1626, HCC-44 817, SJSA1 773, T-47D 697, OVTOKO 627, LN-229 579, YKG1 531, A-673 164; lineages Lung/CNS/Bone/Breast/Ovary). TRPS1 RNA is effectively a T-47D indicator (top marker of RNA cluster 4). So rho(TRPS1 RNA, GATA-motif deviation) across all 5814 cells is little more than "GATA-motif accessibility in the 12% breast cells vs the 88% non-breast majority." The negative sign is fully explained by GATA-motif deviation being relatively higher in the non-breast lines — consistent with the GATA6-positive lung lines (GATA6 rho +0.139), which dominate the pool. The submission concedes it cannot restrict to T-47D; that concession is not a caveat, it is the entire result.

(4) EFFECT SIZE IS AT THE NOISE FLOOR. rho=-0.148 is ~2% of variance, inside the same |rho| 0.04-0.15 band as CTCF (+0.066), ESR1 (+0.082), FOXA1 (+0.108), GATA4 (+0.045). Decisive internal control: GRHL2 — also a T-47D cluster-4 marker (logFC 3.35) in the same pool, same cluster, same design — gives +0.373 (p=4.6e-191), 2.5x larger and positive. So the pooled design does NOT force weak/negative self-motif correlations; a genuinely self-tracking TF produces a much stronger positive one. TRPS1's value does not stand out from noise.

(5) THE CORROBORATION IS A NON-SEQUITUR. I reproduced the marker evidence exactly (top_cluster_markers(4): TRPS1 first by DE score, logFC 4.15, padj 6.9e-275; cross_modal_marker_check(TRPS1): is_rna_marker true, rna_cluster 4, matched_atac_cluster 4, gene_activity_confirms true). It is real, but it only establishes that TRPS1 is breast-cell-type-specific — it carries no information about the sign of a pooled self-motif correlation, so it cannot corroborate the mechanistic claim. Note also the "ranked above ESR1, PGR and GATA3" framing is a ranking-metric artifact: TRPS1's logFC (4.15) is LOWER than ESR1 (5.01), PGR (5.32) and GATA3 (4.73); it leads only on the DE z-score, which rewards low within-cluster variance.

(6) Minor: get_qc_summary reports 10.6% doublets in a multiplexed 8-line pool, which adds cross-line RNA/ATAC mixing that further degrades pooled per-cell correlations. Data quality is otherwise good (FRiP 0.62, TSS 14.6, ARI 0.725), so QC is not the problem — the inference is.
- Already known: True -- Every biological component of the proposal is already published, and the strongest paper reports it in precisely the relevant system with far better methods.

PMID 38377146 (PLoS Genetics 2024), abstract read: opens by calling TRPS1 "the repressive GATA-family transcription factor (TRPS1)," notes "luminal breast cancer cell lines are particularly sensitive to TRPS1 knockout," and reports via an inducible degron in a luminal breast cancer line that "TRPS1 directly regulates chromatin structure," redistributing ER across the genome. That is the proposed finding — a repressive GATA-family factor shaping chromatin/GATA-element accessibility in a luminal breast line — already demonstrated with acute degradation and direct chromatin readouts, which is vastly stronger evidence than a rho of -0.15.

PMID 30563971 (Oncogenesis 2018), abstract read (the submission's own citation): TRPS1 is an atypical GATA factor that "guides the machinery to specific target sites by recognizing GATA elements" and recruits CHD4/NuRD(MTA2) to repress, including enhancer decommissioning at TP63. The submission cites this as background yet the proposed "mechanistic candidate" adds nothing to it.

PMID 38647255 (Am J Surg Pathol 2024, 19,201 tumors) and PMID 39243111 (Diagn Pathol 2024), abstracts read: TRPS1 is a highly sensitive, routinely used breast-cancer IHC marker, "a nuclear protein highly expressed in breast epithelial cells," with TRPS1+GATA3 co-positivity in 47.4-100% of breast cancers. So "TRPS1 is the top-ranked marker of the breast cluster, above GATA3" is standard diagnostic pathology, not a discovery.

What would be genuinely new is a demonstration that TRPS1 occupancy makes GATA-motif accessibility anti-correlate with GATA3 activity WITHIN luminal cells. The submission does not show that, and the available tools cannot: the correlation cannot be restricted to T-47D, and my GATA2/GATA4/GATA6 results argue against it anyway.

### Finding (confidence: medium)

The T-47D cluster's identity in this dataset is defined at least as strongly by a cytokine-receptor/JAK-STAT hormone axis (PRLR, GHR, ERBB4) as by the classic steroid-receptor axis: PRLR has the largest log-fold-change of any marker of the cluster, exceeding both PGR and ESR1, is cross-modally confirmed, and JAK-STAT receptor signaling comes out as a top GO enrichment of the marker set essentially tied with 'response to estrogen'.

**Evidence:** top_cluster_markers(4, n=30): PRLR logFC 6.46, padj 1.7e-273 — the largest logFC in the top-30 set, above PGR (5.32) and ESR1 (5.01); GHR logFC 4.44, padj 1.1e-190; ERBB4 logFC 3.73, padj 1.4e-226. cross_modal_marker_check(PRLR): RNA marker of cluster 4, matched ATAC cluster 4, gene activity confirms (so the receptor locus is also independently more accessible, not just more transcribed). enrich_gene_set on the top-30 markers: 'Positive Regulation Of Receptor Signaling Pathway Via JAK-STAT (GO:0046427)' adj p = 3.8e-03 (GHR, ERBB4, PRLR), ranked alongside 'Response To Estrogen (GO:0043627)' adj p = 3.8e-03 (KRT19, GATA3, ESR1). Honest caveat: prolactin responsiveness of this line is established background knowledge (e.g. PMID 23410749, title-level only), so what is new here is the quantitative ranking (PRLR > PGR > ESR1 by effect size), the independent ATAC gene-activity confirmation, and the module-level co-enrichment with GHR/ERBB4 — not the existence of PRLR expression itself.

**Judger verdict: struck_down**

- Likely artifact: True -- The underlying numbers replicate exactly, but the INTERPRETATION rests on three invalid inferential steps.

(1) logFC is not an "identity-defining strength" metric in this design. check_for_identity_columns shows this is a pool of 8 genotype-confirmed cell lines (T-47D n=697; the other 7 are lung NCI-H838/HCC-44, osteosarcoma SJSA1, ovarian OVTOKO, glioma LN-229/YKG1, Ewing A-673). So logFC for cluster 4 = T-47D vs. a mean over 7 unrelated non-breast lineages, in all of which PRLR, PGR and ESR1 are ~0. With a near-zero denominator the ordering collapses to relative transcript abundance in T-47D plus pseudocount/dropout behaviour, and it is entirely contingent on an arbitrary comparison set (replace the 7 lines with other ER+ breast lines and PRLR's logFC would collapse). Notably the statistic the list is actually RANKED by tells the opposite story: TRPS1 is rank 1, and ESR1 (padj 8.3e-245) and PGR (6.3e-242) both beat ERBB4 (1.4e-226) and GHR (1.1e-190).

(2) The cross-modal confirmation is non-discriminating. I ran cross_modal_marker_check myself on all three: PRLR, ESR1 and PGR are ALL is_rna_marker=true, matched_atac_cluster=4, gene_activity_confirms=true. Citing PRLR's ATAC confirmation as evidence that the cytokine axis rivals the steroid axis is invalid — the steroid receptors pass the identical test.

(3) The GO "tie" is a Benjamini-Hochberg artifact. My own enrich_gene_set rerun on the same 30 genes gives Response To Estrogen nominal p=1.77e-05 and JAK-STAT nominal p=2.29e-05; they share adj p=0.00377 only because BH assigns tied adjusted values to adjacent ranks. Estrogen is nominally the stronger term. The JAK-STAT hit is also a 3-gene annotation tautology (PRLR and GHR are class-I cytokine receptors by definition, so the term is guaranteed once they are in the list), whereas the ER/luminal program is far more broadly represented across the marker set (ESR1, PGR, GREB1, GATA3, TRPS1, AFF3, KRT19, and at n=60 also STC2, CA12, RERG, XBP1, INPP4B — canonical estrogen-induced/luminal genes, background knowledge).

(4) The one independent functional test I could run contradicts the claim. tf_motif_correlation: STAT5A rho=-0.024, p=0.065 (non-significant); STAT5B rho=-0.079, p=1.4e-09 (significantly NEGATIVE); ESR1 rho=+0.082, p=3.9e-10 (significantly positive). There is no positive evidence of STAT5 regulatory activity, while the ER axis is the one with supportive motif evidence. I weight this as supporting rather than decisive, because the effect sizes are tiny and GATA3 is also negative (rho=-0.079, p=2.0e-09) despite being a bona fide luminal TF — the motif tool is noisy in a pooled multi-lineage dataset. Biologically this is also the expected "receptor present but unstimulated" case: PRLR/GHR/ERBB4 are ligand-dependent receptors and standard culture medium supplies no prolactin or GH (background knowledge), so receptor mRNA/accessibility is not evidence of an active JAK-STAT axis.

Things that do survive scrutiny and are NOT artifacts: cluster 4 really is T-47D (genotype-confirmed metadata field, 697 cells — no small-sample problem; ATAC-vs-truth ARI 0.906 per get_qc_summary), QC is good (median 4258 genes/cell, FRiP 0.62, TSS enrichment 14.6), there is no condition axis to confound (check_for_condition_groups: control-only), and I verified by extending to n=60 that PRLR's logFC of 6.46 is indeed the largest in the top 60. So the measurements are sound; it is the ranking-based interpretation that does not survive.
- Already known: True -- (not provided)


### Limitations

- "Faults detected"/"diagnosis matches" use heuristic free-text classifiers on the agent's open-ended answer, not exact ground-truth string matching.
- Known-biology checklist items are generated fresh each run via literature RAG, not scored against a fixed pre-written ground-truth file -- recall/precision against a fixed checklist is a different, complementary evaluation this report doesn't repeat.
- Gene activity is a noisier, indirect accessibility proxy than direct RNA counts; cross-modal disagreement on a real marker is expected, not necessarily a data-quality issue.
- This report reflects a single run on one model; consistency across repeated runs (or across models) is a separate axis this report doesn't cover.
- ATAC downsampling and doublet injection (used for the public dataset) could not be reused as-is here: this dataset has no raw fragments file (downsampling needs one) and RNA `.X` has no raw counts layer (doublet injection needs one) -- a real format constraint, not a design choice. The fault-injection section uses cell-line-label-swap and shuffled-RNA-ATAC-pairing instead.
- This dataset has no drug/treatment condition axis (checked by column name via `check_for_condition_groups`, not assumed) -- every cell is a control. The pipeline supports condition-level analysis (per-arm QC via `condition_group_qc`, condition-specific literature RAG) for a future dataset that does have one; it's simply not exercised here.
- The negative-control check was NOT run for this report. See the public report's own negative-control section for what this check looks like when it does run.


---

## Cost breakdown

- loader_decision: $0.0356
- checklist_generation: $10.5883
- novelty_proposal: $6.2832
- judging: $12.0344
- fault_injection: $13.4132
- **Total: $42.3547**
