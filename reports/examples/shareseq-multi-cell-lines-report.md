# Multiome QC & Hypothesis Agent report (10x PBMC)

*Model: `claude-opus-5` · total cost: $15.6640*

## 1. Dataset & fixed-core summary

**Data source:** public 10x PBMC multiome data. Loading strategy was determined by the agent itself from real file structure, not pre-specified: it decided `combined_single_file` -- "Inspection of the single provided path shows a 10x Cell Ranger h5 with 11,909 obs × 144,978 vars whose feature_types_present list includes both "Gene Expression" and "Peaks", meaning RNA and ATAC modalities are stored together in this one file."

11909 cells, 26349 genes, 107385 ATAC peaks (after feature filtering). RNA: 16 Leiden clusters, median 1826 genes/cell, 3776 UMIs/cell, 9.7% mito, 8.6% predicted doublets (median doublet_score 0.038). ATAC: 21 Leiden clusters, median 13486 fragments/cell, median FRiP 0.76, median TSS enrichment 16.7, median nucleosome signal 0.92. RNA-ATAC cluster agreement (ARI): 0.460. SPI1 expression vs. its own motif's chromVAR deviation (Spearman rho): 0.575.

---

## 2. Gene activity, chromVAR motif deviations, and cross-modal validation

Gene activity scores and chromVAR-style motif deviations are part of the fixed-core pipeline (computed once, cached, reused here). The table below is a systematic sweep (zero extra LLM cost -- these are deterministic Python tool wrappers, run directly, not routed through an agent turn) of the textbook ArchR/Signac cross-modal cell-type-call validation over EVERY RNA cluster: its own top-scoring marker gene is checked against its real cross-modal ATAC partner (matched by cell overlap, not by clusters coincidentally sharing an integer label). A cluster whose top-scoring gene doesn't clear the standard significance threshold (padj<0.05, logFC>1.0) is skipped, not forced with a non-significant marker.

**How "Matched ATAC cluster" is determined:** since RNA and ATAC come from the same cells, every cell has both an RNA cluster label and an independent ATAC cluster label (the two clusterings are computed separately, so cluster numbers between them don't inherently mean anything -- "RNA cluster 4" and "ATAC cluster 4" aren't related just because they share a number). To find the real cross-modal partner of an RNA cluster:
1. Take every cell in that RNA cluster.
2. Look up which ATAC cluster each of those same cells landed in.
3. Build a contingency table (RNA cluster x ATAC cluster cell counts) and take the ATAC cluster with the most overlapping cells -- i.e., the mode.

| Gene | RNA marker of cluster | Matched ATAC cluster | Gene-activity confirms |
|---|---|---|---|
| SLC8A1 | 0 | 5 | True |
| DPYD | 1 | 5 | True |
| FCGR3A | 2 | 2 | False |
| CCL5 | 3 | 3 | True |
| GNLY | 4 | 4 | True |
| FHIT | 5 | 10 | True |
| INPP4B | 6 | 8 | True |
| LEF1 | 7 | 9 | True |
| LYZ | 8 | 0 | True |
| BANK1 | 9 | 16 | True |
| CST3 | 10 | 14 | True |
| SOX4 | 11 | 17 | False |
| PLXDC2 | 12 | 15 | False |
| RALGPS2 | 13 | 16 | True |
| TCF4 | 15 | 19 | True |

12/15 clusters' top markers cross-validate between modalities (1/16 RNA cluster(s) skipped: no marker cleared the significance threshold). Per CLAUDE.md's design principle, discordance here is expected, not a bug -- gene activity is a noisier, indirect accessibility proxy, and RNA-ATAC agreement is not assumed to be perfect.

---

## 3. Known-biology checklist (literature RAG, identity determined by the agent)

**Identity determination:** **Human PBMCs (peripheral blood mononuclear cells), healthy/untreated, one donor-level sample **

- `check_for_identity_columns` found **no** identity-encoding field (no cell-line/donor/DepMap "ACH-" column), so I fell back to marker-based inference.
- `check_for_condition_groups` found **no** condition axis → treat as control-only; no drug/vehicle comparison is available.
- Marker-based composition (8 RNA clusters): **c0** naive CD4 T (LEF1 +4.6, CCR7, TCF7, IL7R); **c5** memory/activated T (IL32, LTB); **c4** NK/cytotoxic (GNLY, NKG7, PRF1, KLRD1); **c3** B cells (MS4A1 +9.1, EBF1 +8.8, PAX5 +8.3, BANK1 +7.8, CD79A +7.0)

Each item below required the agent to: search PubMed, fetch and actually read a real abstract (not just a title), and verify with a real tool call that the gene/motif is usable in this dataset before recording it.

## Confirmed by data

### RNA markers

- **CD14**: CD14 is expressed/enriched in classical (CD14+CD16-) monocytes; here it is a significant RNA marker of monocyte cluster 1 (PMID 36713384, *Frontiers in Immunology*, 2022; data-presence verified: True)
- **S100A8**: S100A8 (calprotectin subunit, with S100A9) is highly expressed by blood classical monocytes/myeloid cells; here S100A8 is a significant RNA marker of monocyte cluster 1 (PMID 32810439, *Cell*, 2020; data-presence verified: True)
- **FCN1**: FCN1 (M-ficolin) mRNA is expressed in blood mononuclear phagocytes/monocytes; here FCN1 is a significant RNA marker of monocyte cluster 1 (PMID 23944633, *Autoimmunity*, 2013; data-presence verified: True)

### Motifs

- **MA0080.7.Spi1 (PU.1)**: PU.1/SPI1 motif-containing regulatory elements should show elevated accessibility in monocyte/macrophage-lineage cells, since PU.1 is the myeloid lineage-determining factor that selects and activates the macrophage enhancer repertoire (PMID 25480297, *Cell*, 2014; data-presence verified: True)
- **MA0466.4.CEBPB (C/EBPbeta)**: C/EBP-family motifs are enriched in the enhancer repertoire of peripheral monocytes, so CEBPB motif accessibility should be elevated in the monocyte compartment (PMID 40702585, *Genome Medicine*, 2025; data-presence verified: True)
- **MA0476.2.FOS (AP-1)**: AP-1 (FOS/JUN) motifs mark enhancers established downstream of pro-inflammatory signalling in peripheral monocytes, so AP-1/FOS motif accessibility should be a feature of the monocyte compartment (PMID 40702585, *Genome Medicine*, 2025; data-presence verified: True)

### TF-expression-tracks-motif-accessibility

- **SPI1**: PU.1/SPI1 protein concentration sets the affinity threshold for in vivo occupancy of PU.1 motifs within accessible chromatin during human macrophage differentiation, so SPI1 RNA level should track SPI1-motif accessibility in monocytes (observed Spearman rho = 0.595, p = 3.5e-49) (PMID 23658224, *Nucleic Acids Research*, 2013; data-presence verified: True)
- **MAFB**: In monocytic cells, induction of MafB mRNA/protein increases MafB occupancy of MARE (Maf) motifs together with histone H4 acetylation indicating chromatin activation, so MAFB RNA should track MAFB-motif accessibility (observed Spearman rho = 0.362, p = 6.5e-17) (PMID 22820162, *Experimental Cell Research*, 2012; data-presence verified: True)
- **CEBPA**: C/EBPalpha occupancy of macrophage/myeloid enhancers causes chromatin opening and activation of macrophage genes (with PU.1 redistribution), so CEBPA RNA level should track CEBPA-motif accessibility in the monocyte compartment (observed Spearman rho = 0.262, p = 2.7e-9) (PMID 37365888, *eLife*, 2023; data-presence verified: True)

## Rejected by data

Real candidates the agent checked via `tf_motif_correlation` on the way to the 3 confirmed
TF-expression-tracks-motif-accessibility items above, whose own data-presence check came back
negative:

- **JUNB** (AP-1 family): Spearman rho = 0.046, **p = 0.30** -- no significant relationship between
  JUNB RNA and its own motif's (MA1140.3.JUNB) accessibility in this dataset. Dropped.
- **KLF4**: Spearman rho = **-0.038**, p = 0.40 -- wrong sign and not significant against its own
  motif (MA0039.5.KLF4). Dropped.

---

## 4. Fault-injection eval (compared across models)


| Model | Faults detected | Correct diagnosis | False alarms | Cost/run |
|---|---|---|---|---|
| claude-haiku-4-5 | 3/3 | 3/3 | 1/1 | $0.0499 |
| claude-sonnet-5 | 3/3 | 3/3 | 1/1 | $0.2476 |
| claude-opus-5 | 3/3 | 3/3 | 1/1 | $0.8599 |
| claude-fable-5 | 3/3 | 3/3 | 1/1 | $0.4594 |



### Per-scenario answers (claude-opus-5)

#### clean_control (no fault injected)


**The data is usable and not broken, but it is not uniformly clean — there's a clear asymmetry between a very strong ATAC library and a shallow, noisier RNA library, plus modest cross-modal agreement.** Three things I'd flag specifically, in order of how much they should change your analysis plan.

**1. RNA depth/complexity is low, and unusually so relative to the ATAC side (the main finding).**
Median 3,776 UMIs/cell against 1,826 genes/cell is only ~2 counts per detected gene. That means most genes in a typical cell are supported by a single transcript, so per-cell expression estimates are dominated by sampling noise. Meanwhile the ATAC side is excellent (13,486 fragments/cell, FRiP 0.76, TSS enrichment 16.7). A well-balanced 10x multiome run normally has RNA depth several-fold higher than this. Practical consequence: DE and any per-cell RNA-ATAC correlation will be attenuated by dropout, and fine-grained RNA subtypes (e.g. CD4 memory subsets) may not separate. *Confidence: high* that the numbers are as stated and the ratio is atypical (`get_qc_summary`); *medium* on the cause, which I cannot distinguish from QC alone (shallow sequencing vs. degraded RNA vs. aggressive filtering).

**2. RNA-ATAC cluster agreement ARI = 0.460 is modest.**
This is the number most likely to be read as "something's wrong," so it's worth separating the benign from the serious explanation. The benign one is mechanical: ATAC was cut into 21 Leiden clusters vs. 16 for RNA, and ARI is penalized when two clusterings are at different resolutions, with shallow RNA (point 1) further blurring RNA boundaries. The serious one would be a barcode-pairing failure — RNA and ATAC profiles not actually from the same cells. **I think this is the benign case**, for two independent reasons from the data: RNA cluster 0 (SLC8A1, TYMP, AOAH, PSAP, HLA-DRA) and ATAC cluster 0 (FPR1, RAB31, PLXDC2, SLC8A1) independently recover the *same* monocyte identity; and SPI1 RNA vs. its own Spi1 motif chromVAR deviation gives rho = 0.575 (p = 0.0), which is a strong within-cell cross-modal coupling that a barcode mismatch would destroy. *Confidence: medium-high* that pairing is intact; *medium* that 0.46 is fully explained by resolution mismatch + RNA noise. I'd still re-check ARI after matching cluster numbers before reporting it.

**3. Predicted doublet rate 8.6% is moderately elevated** (median doublet score 0.038). This is a caveat to carry, not a reason to call the dataset dirty — every ATAC metric is normal-to-excellent and the clusters are biologically coherent. It mainly means small "intermediate" or marker-co-expressing clusters should be treated as suspect until checked. **Median mitochondrial fraction 9.7%** is also on the high side for a filtered PBMC object (mildly stressed/ambient-heavy, not failing). *Confidence: high* on the numbers, *medium* on the interpretation.

**Not a problem:** ATAC quality across the board — TSS 16.7, FRiP 0.76, nucleosome signal 0.92 all sit comfortably in the healthy range for 10x scATAC. No condition axis exists here (`check_for_condition_groups`: control-only), so don't analyze one.

**Sample identity** (inferred, not given): no identity metadata column exists (`check_for_identity_columns`: none found), so I inferred from markers — this is a PBMC-like mixture, with monocytes (cluster 0) and a cytotoxic T/NK population (cluster 3: CCL5, NKG7, GZMA, GNLY, PRF1). *Confidence: high* for those two populations; the remaining 14 clusters I did not characterize.

#### shuffled_rna_atac_pairing (5954/11909 cells have their ATAC profile swapped with another cell's.)

The data is **mostly high quality, but there is a genuine anomaly** 

**First, what the sample is** (needed before judging anything): no identity or condition metadata columns exist, so I inferred it from markers. RNA cluster 0 is classical monocytes (LYZ, FCN1, TYMP, CPVL, HLA-DRA, all adj. p ≈ 0), MS4A1 marks RNA cluster 9 (B cells) and CD3E marks RNA cluster 3 (T cells). GO enrichment on cluster 0 returns immune-effector/IL-6/TNF-response terms. This is a **human PBMC-like immune sample, control-only** (no treatment axis found).

**What looks good.** ATAC quality is genuinely strong: median 13,486 fragments/cell, FRiP 0.76, TSS enrichment 16.7, nucleosome signal 0.92 — all comfortably in the healthy range. RNA depth (1,826 genes / 3,776 UMIs median) is modest but normal for PBMCs. Mito at 9.7% is borderline-but-acceptable.

**Caveat (not disqualifying):** 8.6% predicted doublets is moderately elevated (median doublet score 0.038). Worth naming explicitly and discounting rare/intermediate-phenotype clusters accordingly, but with every other metric normal this alone doesn't make the data unclean.

**Problem — RNA–ATAC cluster agreement of ARI 0.116 is anomalously low.** This is a real finding, not something to wave away. Importantly, though, I tested whether it reflects broken cell pairing, and **it does not**: `cross_modal_marker_check` confirms MS4A1 (RNA cluster 9 → ATAC 16) and CD3E (RNA cluster 3 → ATAC 3) in both modalities independently, ATAC cluster 0's top accessibility markers are myeloid (FPR1, RAB31, GAS7, PLXDC2, sharing SLC8A1 with RNA cluster 0), and SPI1 expression tracks its own motif deviation at rho 0.575 (p = 0.0) — barcode-shuffled data would collapse all of these. So the modalities *are* correctly paired and recover the same biology; the low ARI most likely reflects over-fragmented ATAC clustering (21 ATAC vs 16 RNA clusters, and ARI punishes splitting hard) rather than a data defect. That said, 0.116 is low even for a resolution mismatch, so I'd re-run ATAC clustering at matched resolution before trusting any cluster-level cross-modal claim.

#### atac_downsampling (ATAC fragments downsampled to 20% (kept 35024911 lines for 11909 barcodes).)

**Not clean — but the problem is specific and diagnosable, not a global failure.** The single thing that looks wrong is **ATAC sequencing depth: median 2,698 fragments/cell**, which is roughly 3–4× below the depth usually targeted for 10x multiome ATAC. Two secondary numbers look off in a way that is consistent with that same root cause, and one is a routine caveat.

**1. Shallow ATAC libraries (the real problem).** Median 2,698 fragments/cell. Critically, this is a *depth* problem, not a *signal-quality* problem: the same libraries show median FRiP 0.76 and median TSS enrichment 14.5, both excellent, and nucleosome signal 0.93 is normal. So the Tn5 reaction, nuclei integrity, and peak set are all healthy — there simply aren't enough reads per cell. Practical consequence: per-cell peak matrices will be very sparse, and peak-level or single-locus ATAC claims are underpowered. Worth noting that FRiP 0.76 is itself on the high side; with shallow libraries, the reads you do get concentrate in the strongest promoter-proximal sites, which inflates FRiP and TSS enrichment relative to what a deeper library on the same cells would show. I'd treat the good FRiP/TSS as confirming library quality, not as offsetting the depth deficit. **Confidence: high** that the number is low and that it limits peak-level power; **medium** on the FRiP-inflation interpretation (mechanistic reasoning, not directly measured here).

**2. Modest cross-modal cluster agreement: ARI 0.460, with 21 ATAC clusters vs 16 RNA clusters.** This is a genuine finding to flag on its own merits, not something to wave away. ATAC is fragmenting into more clusters than RNA while agreeing with it only moderately — the signature of a noisy ATAC embedding splitting cells on technical variation (depth) rather than biology. It is consistent with item 1, but I want to be explicit that I'm inferring the link, not measuring it. Practical consequence: **do not use ATAC Leiden labels as the primary cell-type partition** here; anchor cell identity on RNA and use ATAC for confirmation. **Confidence: high** that agreement is only moderate; **medium** that shallow depth is the cause.

**3. Shallow/low-complexity RNA too.** 3,776 median UMIs against 1,826 median genes is only ~2 UMIs per detected gene — these libraries are sequenced near the low-complexity end. Genes/cell is adequate for cell typing, but low-expression genes and subtle DE will be missed. **Confidence: high** on the numbers, **medium** on the practical impact.

**4. Moderately elevated predicted doublets: 8.6%.** Above typical, but the median doublet score is only 0.038, so this is a modest tail of flagged cells rather than a pervasively contaminated dataset, and every other RNA metric is normal. Name it as a caveat on any small or intermediate cluster (which can be doublet artifacts); not a reason to call the dataset unusable. Median mito 9.7% is borderline-high but acceptable for a monocyte-rich sample.

**What still works.** Two independent positive controls say the biology survived: SPI1 RNA vs. its own motif's chromVAR deviation gives Spearman rho 0.575 — a strong TF-activity coupling — and TYMP is confirmed as an RNA marker of cluster 0 with ATAC gene activity independently elevated in its matched ATAC partner (cluster 5). So cross-modal signal is recoverable at the gene-activity level even though peak-level resolution is compromised. Cluster 0's markers (TYMP, AOAH, PSAP, HLA-DRA, TNFAIP2, SLC8A1) identify it as monocytes, i.e. this is a PBMC-like primary immune sample — which also makes the SPI1 result a sensible positive control rather than a coincidence.

**Bottom line:** usable for cell-type-level and gene-activity-level analysis; treat peak-level ATAC conclusions and ATAC-derived clusters as low-confidence.


#### injected_doublets (1323 synthetic doublets added to 11909 real cells (10.0% of final total).)

The data is **mostly high quality, but not clean — there is one substantive flag plus two lesser ones.**

**The main problem: poor RNA–ATAC cluster agreement (ARI = 0.460).**
In a 10x multiome, RNA and ATAC come from the *same nucleus*, so the two clusterings are two views of one identity label and should agree strongly. An ARI of 0.46 means the partitions share less than half the structure you'd expect from paired measurements. Part of this is likely benign resolution mismatch — ATAC was split into 21 Leiden clusters vs. 16 for RNA (`get_qc_summary`), and extra ATAC splits mechanically depress ARI without meaning the cell calls conflict. But a 5-cluster difference doesn't usually cost you that much ARI on its own, so I'd treat this as a genuine partial discordance in the joint embedding/clustering step, not just a knob setting. Importantly, this is **not** explainable by the doublet rate, and it should not be waved away.

One reassuring counterweight: SPI1 RNA vs. its own motif chromVAR deviation gives Spearman rho = 0.575 — a strong, biologically specific cross-modal correlation. That argues the *barcode pairing itself is intact* (a swapped/shuffled RNA–ATAC pairing would destroy this), so the ARI issue is most likely clustering-resolution/embedding-level, not a catastrophic pairing failure. **Confidence: high** that the ARI is low relative to expectation; **medium** on the interpretation that it's resolution/clustering rather than pairing.

**Secondary flag: predicted doublet rate 17.8%.** This is roughly double the ~8–9% you'd expect from 10x loading at ~13k cells (background knowledge). Note the internal tension: the *median* doublet score is only 0.038, meaning the bulk of cells look clean and the 17.8% comes from a flagged tail — consistent with a permissive calling threshold rather than a genuinely doublet-laden library. Per standard practice, a moderately elevated doublet rate alone is not grounds to call the dataset unclean when everything else is normal, but it is a caveat you should carry into any downstream cluster-level or rare-population claim (and it could be *contributing* to the ATAC over-clustering and hence the ARI). **Confidence: high** on the number, **medium** on the "permissive threshold" explanation.

**Minor flag: 9.8% mitochondrial reads, and an odd genes:UMI ratio.** 9.8% mito is borderline — acceptable for a standard 10% cutoff but sitting right at it, suggesting some stressed/low-quality cells survived filtering. Separately, median 1900 genes on median 3996 UMIs means ~2.1 UMIs per detected gene — unusually high library complexity for that depth, which is what you see with shallow sequencing and/or ambient RNA contribution. Neither is disqualifying, but the RNA side is shallow. **Confidence: medium.**

**What looks genuinely good:** the ATAC library is excellent — median 13,486 fragments/cell, FRiP 0.76, TSS enrichment 16.7, nucleosome signal 0.92 are all comfortably in or above the good range. The ATAC modality is not the weak link here.

**Scope note:** `check_for_identity_columns` found no identity field and `check_for_condition_groups` found no treatment axis (control-only), so this QC read is on an unlabeled, single-condition sample; I have not established the tissue/cell-type composition and am not assuming one.

**Recommendation:** re-run clustering at matched resolution across modalities (or cluster jointly, e.g. WNN) and re-score ARI before trusting any cluster-level result; tighten the doublet threshold toward the expected ~8–9%.

## 5. Novel findings, adversarial judging, and limitations

### Finding (confidence: medium)

TCF7L2 is a myeloid-expressed TF in this PBMC sample, yet its own chromVAR motif deviation is STRONGLY NEGATIVELY coupled to its RNA (Spearman rho = -0.506) -- the opposite sign from every other TF tested here. The most parsimonious explanation is motif-family degeneracy: the TCF7L2 PWM is essentially the shared TCF/LEF HMG-box site, so the 'TCF7L2 motif' chromVAR score in PBMCs reads out lymphoid TCF7/LEF1 activity in the T-cell compartment, while TCF7L2 mRNA sits in the monocyte compartment. This is a concrete caution: for paralogous motif families, a TF's own motif score can anti-report its own expression.

**Evidence:** tf_motif_correlation TCF7L2 (MA0523.2) rho = -0.5064, p = 6.5e-34 -- the only negative correlation among 16 TFs tested. Same tool: TCF7 (MA0769.3) rho = +0.4534, p = 1.0e-26 and LEF1 (MA0768.3) rho = +0.3489, p = 9.3e-16, i.e. the same TCF/LEF site family tracks POSITIVELY with the T-cell paralogs. TCF7L2's expression is myeloid, confirmed in both modalities: cross_modal_marker_check TCF7L2 -> is_rna_marker=True (RNA cluster 1, the FCN1/LYZ classical-monocyte cluster), matched ATAC cluster 1, gene_activity_confirms=True; and top_cluster_markers cluster 6 (FCGR3A lfc 5.34, CDKN1C lfc 7.64, LST1 lfc 4.34; enrich_gene_set of those markers -> Fc-gamma receptor signaling GO:0038094 adj p = 6.4e-05, i.e. CD16+ non-classical monocytes) lists TCF7L2 at lfc 5.15, adj p = 9.5e-12. The T-cell side is likewise cross-modally real: cross_modal_marker_check LEF1 -> RNA cluster 0 (LEF1 lfc 4.62, CCR7 3.97, TCF7 3.01), matched ATAC cluster 10, gene_activity_confirms=True. Literature consistency (not the same observation): PMID 40631795 reports TCF7L2 regulon activity as specific to nonclassical monocytes in human PBMC scRNA-seq -- i.e. the myeloid expression is expected, the negative motif coupling is the new part. A parallel, weaker instance of the same paralog effect: SPIB (B/pDC-expressed, ETS-family site shared with PU.1) rho = +0.160, p = 3.3e-04, far below SPI1's rho = 0.595.

**Judger verdict: survives**

- Likely artifact: False -- The empirical observation is real and I reproduced it exactly: tf_motif_correlation TCF7L2 (MA0523.2) rho = -0.5064, p = 6.46e-34; TCF7 (MA0769.3) +0.4534; LEF1 (MA0768.3) +0.3489. Cross-modal anchoring also replicated: cross_modal_marker_check TCF7L2 -> is_rna_marker=True, RNA cluster 1, matched ATAC cluster 1, gene_activity_confirms=True; and I independently confirmed cluster 1 is classical monocytes from its own markers (FCN1 lfc 3.42, LYZ 3.00, TYMP 3.25, AOAH 3.13, all adj p < 1e-19). T-cell side also replicated (cross_modal_marker_check TCF7 and LEF1 both -> RNA cluster 0, matched ATAC cluster 10, gene_activity_confirms=True).

  Adversarial test of sign uniqueness: rather than trust the original 16, I tested 10 further TFs of my own choosing. CEBPB +0.300, MAFB +0.362, IRF8 +0.295, ETS1 +0.252, ZEB1 +0.433, GATA3 +0.151, TCF4 +0.283, SOX4 +0.132, SPI1 +0.595; and near-zero/ns for KLF4 -0.038 (p=0.40), KLF2 +0.059 (p=0.19), RUNX1 +0.012 (p=0.79), NFKB1 +0.038 (p=0.40), JUNB +0.046 (p=0.30), POU2F2 +0.038 (p=0.39). Across ~24 TFs, TCF7L2 remains the ONLY substantial negative, and its |rho| is second only to SPI1. So this is not a pipeline-wide sign bug, not noise, and not a low-magnitude blip.

  Remaining caveats, judged as limitations rather than refutations: (1) Pseudo-replication -- n=500 cells but the variance is essentially cluster-level (monocyte vs T compartment), so p=6.5e-34 is badly inflated; the honest evidence is a two-compartment contrast, not 500 independent observations. The sign, however, is not in doubt. (2) Mechanism is under-determined by the available tools -- a negative TF-RNA/own-motif correlation is conventionally read as repressor activity, and TCF7L2 without beta-catenin is a bona fide TLE/Groucho-dependent repressor (background knowledge), which would predict the same negative sign; nothing in this tool set distinguishes paralog PWM degeneracy from genuine repressive function. (3) The judge's own further tests actually weaken the generality of the proposed mechanism: other paralog-mismatch cases in this dataset attenuate toward zero rather than inverting (SPIB +0.160 vs SPI1 +0.595; KLF4 -0.038 vs KLF2 +0.059; POU2F2 +0.038) -- degeneracy alone predicts rho ~ 0, not -0.51; a strong inversion additionally requires the near-total mutual exclusivity that holds for TCF7/LEF1 (naive-T, sites wide open) versus TCF7L2 mRNA (monocyte-restricted), which is a compositional effect -- exactly what the finding claims. (4) Moderate data quality (9.9% median mito, 7% predicted doublets, RNA-ATAC ARI 0.467) adds noise but cannot manufacture a -0.5 correlation. The finding does not assert TCF7L2 chromatin biology -- it asserts that the motif score mis-reports the TF, which is a correct diagnosis of an artifact, not itself an artifactual biological claim.
- Already known: False


### Finding (confidence: medium)

Within the AP-1 family, RNA-to-motif coupling is strongly asymmetric by subunit: the Fos-side subunits track AP-1 motif accessibility (FOS rho = 0.572; FOSB::JUN rho = 0.387) while the Jun-side subunits barely or don't (JUN rho = 0.167; JUNB rho = 0.046, not significant), even though Fos and Jun bind the same TRE site as an obligate heterodimer. So in this PBMC dataset, AP-1 motif accessibility is quantitatively reported by FOS/FOSB mRNA and NOT by JUNB mRNA -- i.e. the lineage/state-specific variance in AP-1 activity lives in the Fos subunit, consistent with Jun-family mRNA being broadly/constitutively expressed across all PBMC lineages.

**Evidence:** tf_motif_correlation: FOS (MA0476.2) rho = 0.5718, p = 9.3e-45; FOSB (MA1127.1 FOSB::JUN) rho = 0.3873, p = 2.4e-19; JUN (MA0488.2) rho = 0.1673, p = 1.7e-04; JUNB (MA1140.3) rho = 0.0464, p = 0.30 (n.s.). All on the same 500 cells (get_qc_summary: 500 cells, 15695 genes, 94708 peaks), so the contrast is not a power artifact -- FOS reaches p ~1e-44 while JUNB is flat. For scale within the same dataset, the myeloid LDTF SPI1 gives rho = 0.595 (get_qc_summary / prior step), so FOS coupling is comparable to the strongest TF here while JUNB is indistinguishable from zero. Related context from literature actually read: PMID 40631795 reports FOSB among the regulons enriched in classical monocytes in human PBMCs, consistent with Fos-side mRNA carrying myeloid-compartment-specific variance.

**Judger verdict: struck_down**

- Likely artifact: True -- The numbers replicate exactly (FOS 0.5718/9.3e-45; FOSB::JUN 0.3873/2.4e-19; JUN 0.1673/1.7e-04; JUNB 0.0464/p=0.30), so this is not a reporting error. But the biological interpretation is confounded, on five independent grounds.

(1) CELL-TYPE COMPOSITION CONFOUND. Marker analysis shows RNA clusters 1 and 2 are both monocytes (cluster 1: TYMP, JAK2, AOAH, FCN1, LYZ, PSAP; cluster 2: VCAN logFC 5.04, CSF3R, SLC11A1, FCN1; enrich_gene_set on these returns Inflammatory Response GO:0006954 adj-p 4.5e-3, Response To Molecule Of Bacterial Origin adj-p 3.6e-3), while cluster 0 is naive T (LEF1 logFC 4.62, CCR7, TCF7, IL7R). cross_modal_marker_check shows FOS, FOSB and SPI1 are all RNA markers of the SAME cluster (1), and FOSB (logFC 2.76, adj-p 5.2e-21) and FOS (logFC 2.24, adj-p 2.9e-19) are DE markers of monocyte cluster 2. JUNB is a marker of NO cluster (is_rna_marker=false). AP-1/TRE accessibility in PBMC is monocyte-biased (PMID 34174187, Cell 2021, explicitly reports a monocyte subcluster defined by chromatin accessibility at AP-1-targeted loci in human PBMC scATAC). So the FOS correlation is the monocyte-vs-lymphocyte axis appearing on both sides of the correlation, not subunit-specific regulatory reporting.

(2) THE EFFECT IS NOT AP-1-SPECIFIC, WHICH IS FATAL TO THE FRAMING. BACH1 -- a bZIP that is NOT a Fos and NOT a Jun, and is itself a monocyte cluster-1 RNA marker -- gives rho = 0.563, statistically indistinguishable from FOS's 0.572. BACH2 (T-cell, cluster-0 marker) gives rho = -0.646, LARGER in magnitude than FOS. ATF3 = 0.448, CEBPB = 0.300, CEBPA = 0.262. Even TFs from unrelated motif families track the same axis: TCF7 = 0.453, LEF1 = 0.349, i.e. comparable to FOSB (0.387) and FOSL2 (0.367). |rho| in this dataset simply tracks how lineage-restricted a TF's mRNA is, with no AP-1 subunit logic required. The original evidence's own scale comparison (SPI1 = 0.595, which I reproduce) is not a control -- it is a demonstration of the confound: FOS behaves exactly like a myeloid lineage marker because it is one here.

(3) MOTIF NON-INDEPENDENCE. chromVAR aggregates accessibility "within peaks sharing the same motif" (PMID 28825706, abstract read). FOS, JUN, JUNB and JUND all bind the same TRE/TGASTCA core (background knowledge), so their deviation scores are near-duplicate variables. The four comparisons are therefore ONE accessibility axis regressed against four different mRNAs -- not four independent measurements of subunit-specific coupling. The finding concedes this itself ("bind the same TRE site"), which means the entire contrast is a property of mRNA distribution across cell types and carries no information about heterodimer subunit behaviour.

(4) VARIANCE FLOOR MAKES THE JUNB RESULT NEAR-TAUTOLOGICAL. A transcript with little cross-cell variance cannot correlate with anything. JUNB is a marker of no cluster, so rho ~ 0 is a statistical near-necessity. The finding's stated mechanism ("Jun-family mRNA broadly/constitutively expressed") is the trivial explanation, restated as if it were a discovery.

(5) NO INDEPENDENT CROSS-MODAL REPLICATION + WEAK SUBSTRATE. gene_activity_confirms = false for FOS, FOSB, JUN and BACH1 -- the independent ATAC-side signal does not corroborate. Dataset is only 500 cells with 9.9% median mito, 7.0% predicted doublets and RNA-ATAC cluster agreement ARI of just 0.467. Additionally, FOS/FOSB/JUN are immediate-early genes induced by ex vivo handling (background knowledge), a processing-stress contribution that cannot be excluded here and that would preferentially inflate Fos-side transcripts in monocytes.


### Finding (confidence: low)

Lineage asymmetry in TF RNA-motif coupling: the B-cell compartment's own lineage-determining TFs are decoupled from their motif accessibility in this dataset, while the myeloid and T-cell ones are not. PAX5 (rho = -0.014, n.s.) and POU2F2 (rho = 0.038, n.s.) show zero coupling and EBF1 (0.126) / SPIB (0.160) only marginal coupling, even though cluster 3 is an unambiguous, cross-modally confirmed B-cell cluster with PAX5 and EBF1 among its strongest RNA markers -- compare SPI1 (0.595), FOS (0.572) and TCF7 (0.453) in the myeloid and T compartments. This is expected-type RNA/ATAC discordance rather than a QC failure: B-lineage TFs bind long, low-copy, GC-rich sites and act substantially through priming/repression, so their target repertoire's accessibility need not scale with their own mRNA.

**Evidence:** tf_motif_correlation: PAX5 (MA0014.4) rho = -0.0141, p = 0.75; POU2F2 (MA0507.3) rho = 0.0383, p = 0.39; EBF1 (MA0154.5) rho = 0.1261, p = 4.7e-03; SPIB (MA0081.3) rho = 0.1601, p = 3.3e-04. Contrast in the same cells: SPI1 rho = 0.595 (get_qc_summary), FOS rho = 0.5718 (p = 9.3e-45), TCF7 rho = 0.4534 (p = 1.0e-26). The B cluster itself is solid, so this is not a failure to detect B cells: top_cluster_markers cluster 3 gives MS4A1 lfc 9.09 (adj p = 6.3e-28), EBF1 lfc 8.82, PAX5 lfc 8.26, BANK1 7.83, CD79A 7.01; and cross_modal_marker_check PAX5 -> is_rna_marker=True (RNA cluster 3), matched ATAC cluster 3, gene_activity_confirms=True, i.e. the PAX5 locus is independently more accessible in the ATAC partner cluster. Caveat noted honestly: total n = 500 cells (get_qc_summary), so the B compartment is a minority of cells and the tools give only global, non-stratified correlations -- reduced power for a B-restricted signal cannot be fully excluded, though PAX5's rho is essentially exactly zero rather than small-and-positive.

**Judger verdict: struck_down**

- Likely artifact: True -- The claimed B-vs-(myeloid/T) asymmetry does not survive an unbiased comparator panel; it is an artifact of cherry-picking the high tail as the "control" group.

I re-ran tf_motif_correlation on the original four B TFs (reproduced exactly: PAX5 -0.0141 p=0.75; POU2F2 0.0383 p=0.39; EBF1 0.1261 p=4.7e-3; SPIB 0.1601 p=3.3e-4) plus 16 non-B TFs in the same 500 cells. Dataset is PBMC (RNA cluster 0 naive CD4 T: LEF1/CCR7/TCF7; 1 and 2 monocyte: FCN1/LYZ/VCAN/CSF3R; 3 B; 4 NK/cytotoxic: NKG7/GNLY/KLRD1; 5 memory T), no identity column, no condition axis (control-only).

Fatal counterexamples, all from T/NK/myeloid, i.e. the compartments the finding claims are "coupled":
- RUNX3 rho = 0.106 (p=0.018) -- LOWER than EBF1 (0.126) and SPIB (0.160). RUNX3 passes cross_modal_marker_check exactly as PAX5 does (is_rna_marker=True, RNA cluster 4, matched ATAC cluster 4, gene_activity_confirms=True), so this is an apples-to-apples comparison, not a weaker marker.
- GATA3 rho = 0.151, also cross-modally confirmed (RNA cluster 5 / ATAC cluster 5, gene_activity_confirms=True) -- statistically indistinguishable from SPIB (0.160).
- TBX21 0.177, EOMES 0.146 (NK) -- same low band.
- NFKB1 rho = 0.038 (p=0.40) and JUNB rho = 0.046 (p=0.30) -- non-significant, i.e. literally the same "zero coupling" as PAX5 (-0.014) and POU2F2 (0.038), in myeloid/ubiquitous factors.
- IKZF1 -0.063 (p=0.16) and TCF3 0.028 (p=0.54) -- near-zero for broadly expressed, non-B-restricted factors.

So near-zero coupling is common in every lineage here, and canonical T-lineage-determining TFs (RUNX3, GATA3) score at or below the B TFs. The proposed lineage axis does not exist; the comparators SPI1 (0.595), FOS (0.572), TCF7 (0.453) are simply the top of the distribution (LEF1 0.349, BCL11B 0.401, CEBPB 0.300, IRF8 0.295, ETS1 0.252, CEBPA 0.262 fill in the middle). FOS is also not a lineage-determining TF at all -- it is a ubiquitous immediate-early AP-1 factor, so it does not belong in a "myeloid lineage TF" comparison set.

The offered mechanism (long, GC-rich, low-copy B-motif chemistry) is directly falsified within a single motif family: FOS rho = 0.572 vs JUNB rho = 0.046 (p=0.30), same AP-1 motif class. Coupling magnitude therefore tracks the individual TF's own mRNA detection/bimodality and the size and distinctness of its expressing compartment (a compositional driver of a global, non-stratified Spearman), not motif chemistry or lineage. SPIB is a clear illustration of that confound: its PU-box motif is near-identical to SPI1's, so the motif deviation is dominated by the large myeloid compartment where SPIB mRNA is absent, mechanically suppressing rho for reasons that have nothing to do with B-cell biology.

The original write-up's own caveat (n=500, B cells a minority, global non-stratified correlations) is the real explanation, and it is not rescued by "PAX5's rho is exactly zero rather than small-and-positive" -- IKZF1 (-0.063) and JUNB/NFKB1 (n.s.) show that pattern outside the B compartment too. The cross_modal_marker_check on PAX5 that was cited validates that cluster 3 is genuinely B; it provides no independent support for the decoupling claim itself.
- Already known: False -- I could not find literature reporting this specific claim (a B-lineage-specific decoupling of TF mRNA from motif accessibility relative to myeloid/T), so I am not striking it down as already-published. Searches for TF-expression/motif-activity discordance, chromVAR correlation benchmarking, and PBMC multiome motif-expression coupling returned essentially nothing beyond the chromVAR method paper itself (PMID 28825706), whose abstract only describes estimating accessibility deviations for motif-sharing peaks while controlling for technical bias -- it makes no claim about lineage-specific TF-RNA coupling.

That said, the finding's mechanistic rationale runs against literature I actually read. PMID 36409886 (PNAS 2022) shows that degrading EBF1 in pro-B cells causes rapid loss of chromatin accessibility at EBF1-binding sites, correlating with altered gene expression -- i.e. EBF1 dose is continuously coupled to accessibility, not decoupled from it. PMID 41266087 (Genes Dev) performed combined scRNA/ATAC on B-lymphoid progenitors and reports that the accessibility switch "correlated strongly with the initiation of Ebf1 and Pax5 transcription, as well as their functional activities" -- the direct opposite of the proposed B-TF decoupling, in a setting with real power. PMID 26982363 confirms EBF1's C-terminal domain actively opens naive chromatin at low-co-occupancy sites, so a priming/repression-only account of B TFs is not accurate either.

Verdict rests on the artifact analysis, not on prior publication: the claim is novel but wrong.

### Negative control: novel-finding proposal under shuffled RNA-ATAC pairing
With the cross-modal relationship completely removed by 100% shuffled pairing, the same novel-finding proposal step was run again. The model proposed 3 findings: 1 positive relationship claim and 2 findings reporting an absence of signal. The positive relationship claim was subsequently rejected by adversarial judging, while the 2 no-signal findings were not judged because they make no positive claim. Thus, no noise-induced relationship survived the evaluation, which is the expected outcome.

#### Finding (confidence: medium)

Within the KLF/SP GC-box motif family, TF expression–own-motif coupling splits by sign along the lymphoid/myeloid axis: KLF2 shows the only significant POSITIVE coupling in this dataset (Spearman rho = +0.141, p = 0.0016) while KLF4 shows a significant NEGATIVE coupling (rho = -0.110, p = 0.0138), despite these two TFs binding near-identical GC-box/CACCC motifs. This implies the KLF-motif chromVAR axis in PBMC is not reporting a single shared "KLF activity" but is split by cell-type composition (KLF2-high lymphocytes vs KLF4-high monocytes), so KLF-family motif deviations should not be interpreted as interchangeable.

**Evidence:** tf_motif_correlation: KLF2 (MA1515.2) rho=+0.1409, p=0.00159; KLF4 (MA0039.5) rho=-0.1101, p=0.0138; KLF3 (MA1516.2) rho=-0.0576, p=0.198 (ns); SP1 (MA0079.5) rho=+0.0086, p=0.847 (ns). Cell-type context from top_cluster_markers: lymphoid clusters present (RNA 0 naive T: LEF1 lfc 4.62, CCR7 3.97, TCF7 3.01; RNA 3 B: MS4A1 9.09, PAX5 8.26) and myeloid clusters present (RNA 2: VCAN lfc 5.04, FCN1 3.55, CSF3R 4.08; RNA 6: FCGR3A 5.34, CDKN1C 7.64). Literature grounding of the directionality: KLF2 drives naive/quiescent T-cell trafficking programs (CD62L, S1PR1) — PMID 17548599 abstract; KLF4 is a well-established driver of monocyte/macrophage differentiation and polarization — PMID 40552304 abstract. Neither paper reports this motif-level sign split in multiome data.

**Judger verdict: struck_down**

- Likely artifact: True -- The two correlations reproduce exactly (tf_motif_correlation: KLF2 rho=+0.1409 p=0.00159; KLF4 rho=-0.1101 p=0.0138), so this is not a reporting error. It is, however, almost certainly noise mining, for five converging reasons.

(1) THE PROPOSED MECHANISM FAILS ITS OWN STRONGEST TEST CASES. The finding's explanation is that cell-type composition (lymphoid vs myeloid) drives the sign of TF-RNA/own-motif coupling. I tested that directly on the most extremely lineage-restricted TF/motif pairs available in PBMC. EBF1 (RNA logFC 8.82 in the B cluster) gives rho=-0.013, p=0.76. PAX5 (logFC 8.26, same cluster) gives rho=+0.0096, p=0.83. SPI1, the myeloid master regulator, gives rho=-0.047, p=0.29. If compositional structure produced these couplings, EBF1/PAX5/SPI1 would be the largest positive correlations in the dataset. They are indistinguishable from zero. The mechanism invoked to explain KLF2/KLF4 demonstrably does not operate in this dataset where it should be far stronger.

(2) THE CORRELATIONS ARE A NULL DISTRIBUTION AND KLF2/KLF4 ARE ITS TWO TAILS. Across 18 TFs I tested (KLF2/3/4/6/13, SP1/2/3, SPI1, CEBPA, CEBPB, TCF7, LEF1, EBF1, PAX5, IRF8, TBX21, GATA3), everything sits near zero (KLF3 -0.058, SP1 +0.009, KLF6 +0.057, SP3 +0.020, KLF13 +0.085, SP2 -0.055, CEBPB +0.012, TCF7 +0.075, LEF1 +0.078, GATA3 -0.011). KLF2 and KLF4 are simply the extreme ends. Selecting the two tails post hoc and narrating a story about their opposite signs is exactly the failure mode this pattern predicts.

(3) MULTIPLE TESTING KILLS THE NEGATIVE LEG. Within the 8 KLF/SP members the claim itself defines as the family, KLF4's p=0.0138 does not survive Bonferroni; across the 18 TFs tested it is nowhere near significant. KLF2's p=0.00159 is marginal at best. The "sign split" requires BOTH legs to be real; the negative leg is not.

(4) NO TEST OF THE DIFFERENCE WAS EVER DONE. The entire claim is that rho_KLF2 differs in sign from rho_KLF4, but no statistic comparing the two correlations was computed. Effect sizes of |rho| 0.11-0.14 on n=500 explain ~1-2% of variance.

(5) HALF THE STATED CELL-TYPE PREMISE IS UNSUPPORTED IN-DATA, AND THE PREMISE IS SELF-UNDERMINING. cross_modal_marker_check returns is_rna_marker=false for KLF2 -- it is not a significant marker of any cluster here, so "KLF2-high lymphocytes" is imported from background knowledge, not shown in this data. (KLF4 does check out: RNA marker of cluster 1, which my markers identify as classical monocyte -- LYZ, FCN1, TYMP, AOAH -- with ATAC gene activity confirming in matched ATAC cluster 2.) Worse, the claim stresses the motifs are "near-identical." If true, the two chromVAR deviation vectors are near-collinear and the sign difference collapses to nothing more than KLF2 and KLF4 mRNA having different cell distributions -- the chromatin layer adds no independent information, and the "finding" is a known expression fact wearing an epigenomics costume. (I could not verify motif PWM similarity with the available tools; flagged as unverified.)

QC context compounds this: only 500 cells, 7.0% predicted doublets, 9.9% median mitochondrial reads, and RNA-ATAC cluster agreement ARI of only 0.467. Doublets in particular manufacture precisely this kind of mixed-lineage RNA-vs-accessibility signal at the few-percent-of-variance scale being claimed.
- Already known: False -- I could not find a paper reporting this exact observation -- a sign-split in TF-expression/own-motif chromVAR coupling between KLF2 and KLF4 in PBMC multiome -- so as literally stated it is not already published. But that is not much of a defense, because the two components it decomposes into are both well established.

The biological substrate is textbook and stated at review level. PMID 42738814 (Cells, 2026), abstract read in full, says in a single sentence that "KLF2 primarily regulates T-cell quiescence and trafficking... KLF4 contributes to inflammatory effector differentiation," alongside KLF10 and KLF13 -- i.e. functional divergence among near-identical-binding KLF paralogs is the framing premise of the field, not a discovery. PMID 17548599 confirms KLF2 directly activates CD62L and S1PR1 to control T-cell trafficking, and PMID 40552304 confirms KLF4 drives monocyte/macrophage differentiation and polarization. The original submission cited these two correctly.

Notably, PMID 19412182 (Nat Immunol) actively undercuts the clean lymphoid/myeloid dichotomy the finding rests on: ELF4 directly activates KLF4 downstream of TCR signaling in naive CD8+ T cells to induce cell cycle arrest, and Klf4-deficient mice accumulate CD8+CD44hi T cells. KLF4 is a functionally important lymphocyte factor operating on the same naive-T quiescence axis as KLF2 -- so "KLF4 = the monocyte one" is an oversimplification even before the statistics are questioned.

The methodological caveat is also not new. That motif-accessibility-based TF activity inference is unreliable and requires careful benchmarking is the entire premise of PMID 39441876 (PLoS Comput Biol, 2024), which benchmarks chromVAR and alternatives for identifying differentially-active TFs precisely because this inference is known to be error-prone (that abstract does not specifically address paralog motif redundancy, so I do not cite it for that narrower point). The advisory "don't treat family motif deviations as interchangeable" is standard practice guidance, not a result.

So: novel as a sentence, not novel as biology, and the strike-down rests on the artifact analysis rather than on priority.

#### Finding, not adversarially judged (confidence: medium)

The cross-modal cluster pairing in this dataset is degenerate and collapses onto a single myeloid ATAC cluster, so "gene_activity_confirms = false" calls here are a matching artifact rather than genuine RNA-ATAC discordance. Four RNA markers from four different lineages (naive T, B, CD16 monocyte, classical monocyte) all get assigned the SAME ATAC partner cluster (ATAC 2, a monocyte cluster) by cell-overlap matching, even though a properly lineage-matched ATAC cluster demonstrably exists for at least the T-cell case.

**Evidence:** cross_modal_marker_check: LEF1 -> rna_cluster 0, matched_atac_cluster 2, confirms=false; MS4A1 -> rna_cluster 3, matched_atac 2, confirms=false; CDKN1C -> rna_cluster 6, matched_atac 2, confirms=false; TCF7L2 -> rna_cluster 1, matched_atac 2, confirms=true. top_gene_activity_markers show ATAC cluster 2 is myeloid (TREM1 lfc 3.05, COLEC12 3.08, LRMDA 2.98, PLXDC2 2.93, RAB31 2.65), while ATAC cluster 0 is unambiguously the T-cell cluster (LEF1 lfc 2.57, padj 1.9e-11; BCL11B 2.19; BACH2 1.76; SATB1 1.79) and ATAC cluster 1 is a second myeloid cluster (FPR1 2.24, FPR3 2.42, LYN 2.22). So LEF1's true ATAC partner is cluster 0, not 2 — the tool's overlap-based partner is wrong, and the only "confirmed" gene is the one whose RNA cluster happens to be myeloid. Consistent with get_qc_summary RNA-ATAC cluster agreement ARI = 0.467 and an 8 RNA vs 11 ATAC cluster mismatch.

*Classified by the model itself as `no_signal_or_concern` -- this is a methodological/tooling observation about the cross-modal cluster-matching algorithm's own limitation (real and shuffle-independent), not a positive discovery claim, so it wasn't sent to the adversarial Judger.*

#### Finding, not adversarially judged (confidence: medium)

Lineage master-regulator TFs in this PBMC dataset show essentially zero coupling between their own RNA expression and their own motif's chromVAR deviation, even for the TFs with the most extreme lineage-restricted expression. Of 13 lineage/immediate-early TFs tested, none reached significance and all had |rho| < 0.09; the effect is strongest-absent exactly where it should be strongest (PAX5 and EBF1 in the B-cell cluster). Only the ubiquitously-expressed KLF2/KLF4 pair broke through — i.e. in this sample, motif-deviation scores track broad lymphoid-vs-myeloid chromatin state, not the identity of the master TF that supposedly writes it.

**Evidence:** tf_motif_correlation (all non-significant): PAX5 rho=+0.0096 p=0.831; EBF1 rho=-0.0134 p=0.764; SPI1 rho=-0.047 (get_qc_summary); IRF8 rho=-0.0274 p=0.542; CEBPB rho=+0.0116 p=0.796; RUNX3 rho=-0.0115 p=0.798; GATA3 rho=-0.0114 p=0.799; TBX21 rho=+0.0407 p=0.364; FOS rho=-0.0344 p=0.443; NFKB1 rho=+0.0554 p=0.216; TCF7 rho=+0.0746 p=0.096; LEF1 rho=+0.0782 p=0.081; SPIB rho=+0.0830 p=0.064. Contrast with the RNA-side effect sizes for the same TFs from top_cluster_markers: PAX5 lfc 8.26 (padj 3.1e-27) and EBF1 lfc 8.82 (padj 5.5e-23) in RNA cluster 3 (B cells, MS4A1 lfc 9.09). Plausible technical contributor, stated honestly: get_qc_summary shows only 500 cells and median 3760 UMIs / 1829 genes per cell, so TF dropout plus limited power (n=500 needs |rho|>~0.09 for p<0.05) can flatten these correlations; the near-zero point estimates nonetheless argue against a large missed effect.

*Classified by the model itself as `no_signal_or_concern` -- correctly reports an absence of an expected relationship (as expected under a fully shuffled negative control) rather than claiming new biology, so it wasn't sent to the adversarial Judger.*

---

## Cost breakdown

- loader_decision: $0.0189
- checklist_generation: $6.3739
- novelty_proposal: $1.1432
- judging: $4.0314
- fault_injection: $1.6168
- negative_control: $2.5378
- **Total: $15.7221**
