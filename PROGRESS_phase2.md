# Phase 2: peak-to-gene links

Continuation of PROGRESS.md, tracking one analysis-menu addition picked from README's
"Roadmap: future analysis-menu additions" list:

1. **Peak-to-gene links** (cis-co-accessibility, Cicero/ArchR-style) -- a real test for
   the "distal-enhancer regulation" explanation CLAUDE.md's principle 2 has always allowed
   the agent to reach for, but never actually tested. Tests whether a DISTAL peak's
   accessibility tracks the TARGET GENE's own expression (the gene near that peak) --
   never a transcription factor's own expression, which is a different question this
   addition explicitly does not touch.

**Regulon inference was considered and explicitly dropped** (owner correction: regulon
inference as designed still routed through "TF's own RNA level" as one side of the
correlation, which is the exact narrow pattern this phase was supposed to move away from --
see README's existing complaint about novelty proposals collapsing into TF-vs-own-motif).
Cancelled before any regulon code was written; nothing to revert. Peak-level differential
accessibility between two named groups (the README's other roadmap item) was never in scope
for this phase either -- left in the roadmap list.

**No agent/LLM calls were run as part of this phase** -- pure deterministic Python + unit
tests only. The new analysis is wired into the existing tool list so a *future* real run can
use it; this phase doesn't spend any API budget verifying that.

---

## Step 1: Plan

Classified as a **bounded** change (brainstorming skill) -- a new validated analysis wrapper
following the exact conventions already established by `tf_motif_correlation`/
`gene_activity.py`, wired into the existing checklist/novelty agent loops. Explored `core/`,
`menu/`, `agent/loop.py`, `agent/checklist_generator.py`, `agent/novelty.py`,
`agent/prompts.py`, and `tests/` conventions first.

Design approved by the owner, then narrowed (regulon inference dropped -- see note above).
Implementation follows below, one small testable step at a time.

## Step 2: Core analysis -- `core/peak_to_gene_links.py`

Small refactor first: `core/gene_activity.py`'s `_load_protein_coding_gene_windows` used to do
its own GENCODE parsing AND promoter-extension in one function. Split into a public
`load_protein_coding_gene_coords()` (raw, un-extended chrom/start/end/strand per protein-coding
gene symbol) and a thin `_load_protein_coding_gene_windows(upstream_bp)` wrapper that promoter-
extends it -- so peak-to-gene's "proximal exclusion zone" and gene-activity's own window are
built from the exact same parsed gene-body coordinates, not two independent GENCODE-parsing
implementations that could silently drift apart.

`core/peak_to_gene_links.py`: for a gene, finds peaks within `window_bp` (default 500kb, the
standard Cicero/ArchR search radius) of the gene body on the same chromosome, EXCLUDING anything
inside the same promoter-extended window gene activity already counts (so a "link" means
something gene activity doesn't already tell you). For each surviving candidate, Spearman-
correlates its per-cell accessibility (`atac.layers["counts"]`, falling back to `.X`) against the
gene's own per-cell RNA expression -- explicitly never a transcription factor's expression
standing in for the gene, per the owner's correction. BH-corrects (`scipy.stats.
false_discovery_control`, available in the installed scipy 1.17 -- no new dependency needed)
across however many candidates were tested for that gene. Returns the top-N by |rho|, whether any
cleared significance, and the best significant link if any -- a "no link found" result is still a
real, reportable answer (`any_significant_distal_link=False`, `best_link=None`), not an error.

The overlap/exclusion logic (`_select_distal_candidate_peaks`) and the correlation/BH-correction
logic (`_correlate_candidate_peaks`) are both pure functions of plain arrays/frames -- no mdata,
no genome file -- so they're unit-tested on small synthetic coordinate tables
(`tests/test_peak_to_gene_links.py`), same convention as `gene_activity.
_sum_peaks_into_gene_windows`'s existing synthetic tests. Caught and fixed one real bug this way:
an early draft fetched `mdata.mod["atac"]` before checking whether the gene even exists in
`rna.var_names`, so an unknown-gene call crashed with `KeyError: 'atac'` on a bare mock instead of
returning the intended error dict -- moved the `atac` fetch after the existence check.

`menu/peak_to_gene_links.py`: thin wrapper, identical error-dict convention to
`menu/tf_motif_correlation.py`.

**Second real bug, caught by deliberately checking scipy's behavior on degenerate input rather
than trusting the real-data test alone:** `spearmanr` on a zero-accessibility (constant) peak
returns `nan` for both rho and pvalue, and `scipy.stats.false_discovery_control` raises
`ValueError: ps must include only numbers between 0 and 1` on any NaN input -- confirmed directly
(`false_discovery_control([0.01, np.nan, 0.2])` raises). The real end-to-end CD14 test happened
not to hit a zero-variance distal peak, so this would have shipped silently broken for the next
gene that does (sparse, low-count distal peaks are common in ATAC data). Fixed by skipping
zero-variance peaks before correlating (`np.ptp(col) == 0`) and reporting the skip count via a
new `n_skipped_degenerate_peaks` field, rather than crashing the whole gene's lookup. Added
`test_correlate_candidate_peaks_skips_zero_variance_peaks_without_crashing` to lock this in.

**Verified against real data, not just synthetic tests:** `agent_fixed_core_mdata` (the real,
cached small PBMC fixture) was actually present in this environment (the `multiome_agent` conda
env + real raw PBMC data, confirmed via `conda env list` / `data/raw/` contents), so
`test_peak_to_gene_links_real_data_end_to_end` runs the full pipeline for real on CD14 -- real
GENCODE parsing, real proximal exclusion, real per-cell correlation + BH correction -- not just a
mocked/skipped path. All 8 tests in `tests/test_peak_to_gene_links.py` pass
(`/local_home/gual/miniconda3/envs/multiome_agent/bin/python -m pytest tests/test_peak_to_gene_links.py -v`).

## Step 3: Wiring into the agent loop + checklist + novelty

- `agent/loop.py`: added `peak_to_gene_links` to `TOOLS` (right after `tf_motif_correlation`) and
  its dispatch branch in `_execute_tool`. Since `CHECKLIST_TOOLS`/`NOVELTY_TOOLS`/`JUDGER_TOOLS`
  are all `TOOLS + [...]`, the checklist generator, novelty proposer, and Judger all get this tool
  automatically -- no separate wiring needed in those three files beyond the prompt-text changes
  below.
- `agent/prompts.py` (`CORE_SYSTEM_PROMPT`): documented as tool #2 ("twelve tools" total now), and
  principle 2 (discordance) now explicitly tells the agent to call `peak_to_gene_links` rather
  than just naming "distal-enhancer regulation" as a possibility.
- `agent/checklist_generator.py`: added a 4th `CHECKLIST_CATEGORIES` entry, `"peak_to_gene"` -- a
  literature-reported distal regulatory element regulating a named gene, verified via
  `peak_to_gene_links` returning `any_significant_distal_link=true` for that gene in this dataset.
  Updated `CHECKLIST_QUESTION`'s category list and verification-tool mentions accordingly
  ("three categories" -> "four categories" throughout).
- `agent/novelty.py`: added one sentence to `_novelty_question` naming `peak_to_gene_links`
  explicitly as a way to avoid every candidate collapsing into the same TF-vs-own-motif shape --
  directly addressing the README's stated complaint, without touching the regulon idea that was
  dropped.
- `agent/report_generator.py`: added the `"peak_to_gene"` -> `"Distal peak-to-gene links"` label
  and `by_category` entry in `_render_section3`, so a recorded item in the new category actually
  renders instead of silently falling through `labels.get(cat, cat)`'s raw-key fallback.

Ran `pytest --collect-only` across the full suite (202 tests collected, zero import errors) to
confirm every edited file is still syntactically valid and importable, then ran the fast-path
unit tests directly affected by these changes (`test_checklist_generator.py`, `test_novelty.py`,
`test_report_generator.py`, `test_peak_to_gene_links.py`, `test_menu.py` -- 43 passed, all existing
behavior preserved) without waiting on the full real-data suite (`pytest tests/ -q`, which was
also kicked off in the background against the real ~3000-cell fixed-core fixture and real
shareseq-equivalent pipeline -- these are pre-existing, expensive tests unrelated to verifying
this specific change, and take on the order of tens of minutes per PROGRESS.md's own account of
the fixed-core pipeline's cost).

## Step 4: Docs

- `CLAUDE.md`: principle 3's peak-to-gene-links bullet now says "implemented" with a pointer to
  the module/tool, plus a note that the regulon idea was considered and dropped and why. Analysis
  menu bullet now distinguishes implemented items (gene activity / motif-level / peak-to-gene
  links) from the still-unbuilt peak-level differential-accessibility item. Known-biology
  checklist section now lists the 4th category.
- `README.md`: roadmap section rewritten -- peak-to-gene links moved out of "future" into a
  "now built" paragraph with the same explanation; the regulon-inference bullet replaced with a
  short note that it was considered and dropped, pointing here for the full discussion; the
  peak-level differential-accessibility item is unchanged, still the one real gap. "What it
  produces" section 3 blurb now mentions the 4th checklist category by name.

## Step 5: Real-data review surfaced a real statistical gap -- added an effect-size floor

After Step 3, actually ran `peak_to_gene_links` (not just the gated test) on 4 canonical PBMC
markers against the REAL 11,909-cell fixed-core cache (not the small demo cache), to look at what
it found rather than trust the test's structural assertions alone:

| Gene | Candidates | "significant" by q<0.05 ONLY | best \|rho\| |
|---|---|---|---|
| CD14 | 78 | 55 | 0.235 |
| MS4A1 | 72 | 43 | 0.375 |
| CD3E | 127 | 91 | 0.275 |
| NKG7 | 132 | 72 | 0.152 |

This exposed a real design gap: at ~12k cells, p/q-values underflow to ~1e-150 or below for even
a weak rho (~0.1-0.25), so a pure `q<0.05` gate calls most of a gene's genomic neighborhood
"linked" -- not a useful notion of "real distal regulatory element." Fixed by adding a second,
required gate: `significant` now requires BOTH `qvalue < max_padj` AND `|rho| >= min_abs_rho`
(owner-specified threshold: 0.2, after discussing 0.25 first). Re-running the same 4 genes with
the combined gate:

| Gene | Candidates | "significant" (q<0.05 AND \|rho\|>=0.2) | best \|rho\| |
|---|---|---|---|
| CD14 | 78 | 6 | 0.235 |
| MS4A1 | 72 | 5 | 0.375 |
| CD3E | 127 | 1 | 0.275 |
| NKG7 | 132 | **0** | 0.152 |

NKG7 now correctly reports NO significant distal link (`any_significant_distal_link=False`) --
its strongest candidate (rho=0.152) doesn't clear the effect-size floor even though its q-value
alone would have. CD14/MS4A1/CD3E still each find at least one real link. This is a materially
different, more honest result than Step 2 originally shipped with.

Added `DEFAULT_MIN_ABS_RHO = 0.2` to `core/peak_to_gene_links.py`, threaded `min_abs_rho` through
`peak_to_gene_links()` and `_correlate_candidate_peaks()`, and added
`test_correlate_candidate_peaks_effect_size_floor_overrides_tiny_pvalue` -- a synthetic peak with
a real, reproducible rho≈0.16 (found by direct search over `RandomState(1)`/noise-scale
combinations, not guessed) whose q-value alone clears 0.05 but whose effect size doesn't, so the
floor is exercised on a genuine case rather than an edge case that happens to never trigger.
Updated the existing `_correlate_candidate_peaks` test calls to pass the new required argument.
Updated the tool description in `agent/loop.py` and `agent/prompts.py` to state both conditions
explicitly, so the agent (and anyone reading the prompt) knows q-value alone isn't the bar.

All 9 tests in `tests/test_peak_to_gene_links.py` pass, as do the 46 tests across
`test_checklist_generator.py`/`test_novelty.py`/`test_report_generator.py`/
`test_peak_to_gene_links.py`/`test_menu.py` together (the full real-data suite from Step 3 was
still running in the background at time of writing -- pre-existing, expensive, unrelated to this
change).

## Step 6: Answering "what did peak_to_gene_links actually find" with real literature, for free

Asked to check the real significant links found on PBMC (CD14/MS4A1/CD3E/NKG7, see Step 5's table)
against actual literature. Ran `search_pubmed`/`fetch_pubmed_abstracts` DIRECTLY as plain Python
calls (not through `run_agent` / any LLM) -- these hit NCBI's free E-utilities, not the Anthropic
API, so this cost nothing against the owner's API budget, consistent with "don't use my API cost."

Pulled the real genomic coordinates of each significant link and cross-referenced against
`load_protein_coding_gene_coords()` to see what else sits nearby (a peak "linked" to gene X can
really be gene Y's own regulatory activity if it sits inside Y's gene body -- exactly what this
surfaced for CD14). Result: CD3E's link (chr11:118,337,765-118,347,252, 21.6kb away) sits in the
real intergenic gap between CD3D and CD3G (all three in a 50kb cluster on chr11q23) -- fetched and
actually read PMID 2137910 (*Nucleic Acids Research*, 1990), which directly describes a T-cell-
specific "delta promoter and its 3' enhancer; and the epsilon promoter and its 3' enhancer" at
exactly this location. A real, literature-confirmed hit. CD14 and MS4A1's links: several PubMed
searches found nothing specifically describing those exact coordinates (one promising-looking
title turned out, on reading the abstract, to be about a different gene, MD-2 -- caught by
actually reading it rather than citing the title). Also surfaced a real methodological gap: three
of CD14's links sit inside CYSTM1/PFDN1's own gene bodies 300kb+ away -- plausibly those genes' own
open chromatin, not CD14-specific regulation, which `peak_to_gene_links` had no way to flag at the
time. Led directly into Step 7 below.

## Step 7: Added a confound flag -- does a candidate peak sit in some OTHER gene's own window?

Per the owner's follow-up: added `overlapping_genes` to each link in `peak_to_gene_links`'s output
-- for each candidate peak, which OTHER genes' own promoter-extended windows it overlaps (same
window `core.gene_activity` uses for gene activity). Purely informational/confound-flagging; it
does NOT change what gets correlated (still always the queried gene's own RNA, never a TF's,
never another gene's).

Small refactor first to avoid a real perf cost: computing this needed the promoter-extended window
for EVERY OTHER gene in the genome, and the straightforward way to get that
(`core.gene_activity._load_protein_coding_gene_windows()`) re-parses the whole GENCODE annotation
gzip from scratch -- which `peak_to_gene_links` ALREADY does once per call via
`load_protein_coding_gene_coords()` for the queried gene's own coordinates. Calling both would
parse the same multi-MB file twice per call. Split `_load_protein_coding_gene_windows` into a new
public `promoter_extend_gene_windows(coords, upstream_bp)` (pure arithmetic over an
already-loaded coords frame) + a thin wrapper that calls `load_protein_coding_gene_coords()` once
then extends it -- so `peak_to_gene_links.py` can reuse the SAME already-loaded `gene_coords` for
both its own gene's raw coordinates and every other gene's extended window, one parse total.

Added `_find_overlapping_gene_windows` (pure function, chrom-grouped vectorized overlap, same
style as `_select_distal_candidate_peaks`) and 3 new synthetic unit tests, plus wired the result
into `peak_to_gene_links`'s per-link output and the existing real-data end-to-end test. All 9
tests pass.

## Step 8: Regulon inference, take two -- corrected design after a real misunderstanding

The owner clarified (after Step 7) that regulon inference was NOT meant to be dropped for good --
the item being rejected in the original brainstorming round was specifically an AGGREGATE design
(TF RNA vs. one averaged "regulon activity" score across many predicted target genes), which still
reduces to "the TF's own RNA vs. one number" and isn't a claim a literature-checking step can look
up by name. What's actually wanted: a SCENIC+-style regulon analysis whose output is PER-TARGET-GENE
-- for each candidate target, the TF's RNA correlated against that ONE other gene's own RNA, so a
later literature-checking agent gets specific, named claims ("does SPI1 regulate CD14?") instead
of an abstract statistic. Confirmed design scope via `AskUserQuestion` before building (given the
prior back-and-forth on this exact feature): (1) per-target-gene list only, not an aggregate score;
(2) candidate target genes defined by EITHER the TF's motif in the target's own promoter OR a
genuinely significant distal peak-to-gene-style link that itself carries the TF's motif (broader,
more compute, chosen deliberately over promoter-only).

**`core/motif_deviations.py`:** `pc.match_motif` already computes a peak x motif binary match
matrix (`adata.varm['motif_match']`) as a side effect of computing chromVAR deviations, which
`compute_motif_deviations` previously discarded entirely. Added `_persist_motif_match`, which
scatters this matrix (computed on the 2bit-filtered COPY, which can have FEWER peaks than the
real dataset) back onto the FULL peak set as `mdata.mod["atac"].varm["motif_match"]` (scipy
sparse, natively h5-writable, same convention as `gene_activity`'s `obsm` sparse matrices) plus
`mdata.mod["atac"].uns["motif_match_names"]`. Also extracted `best_motif_match(columns, query)` out
of `tf_expression_motif_correlation`'s inline matching logic, so `core.regulon_inference` resolves
a TF's motif the exact same way (exact-token-preferred ranking) against a DIFFERENT name list
without duplicating that logic.

**This means any EXISTING cached fixed-core result (built before this change) lacks `motif_match`**
-- including the agent's own small cache (`agent_fixed_core_n11909_seed0.h5mu`). Rebuilt it with
`force_reload=True` (pure local compute, no API cost -- the same ~agent-fixed-core pipeline this
project already runs, just rerun once) to get real data to test and report against; see Step 9.

**`core/regulon_inference.py`:** two-pass algorithm --
1. *Candidate identification* (no correlation needed to qualify): every peak carrying the TF's
   motif (`atac.varm["motif_match"]` lookup, no re-scanning) that falls in some OTHER gene's own
   promoter-extended window is direct "promoter motif evidence" for that gene
   (`_overlap_pairs`, a pure chrom-grouped vectorized helper). Separately, every TF-motif peak
   within `window_bp` of some OTHER gene (`_genes_near_peaks`, same style, but computing real
   genomic distance rather than a plain overlap test) that is NOT already covered by promoter
   evidence gets its OWN Spearman(peak accessibility, that gene's RNA) test, BH-corrected across
   all such distal pairs for this TF -- only pairs clearing the same `q<0.05 AND |rho|>=0.2` bar
   `peak_to_gene_links` uses become "distal evidence." Genome-wide proximity to a TF-motif peak
   alone is never sufficient.
2. *Reported statistic*: for the UNION of genes qualified either way, Spearman(TF's own RNA, that
   target gene's own RNA), BH-corrected across the candidate set, same significance bar. This --
   not anything from step 1 -- is what `significant`/`spearman_rho` on each target reflects;
   step 1 only decides which genes get tested at all.

**`menu/regulon_inference.py`:** thin wrapper, same error-dict convention as the other two menu
tools.

**Wiring:** added to `agent/loop.py`'s `TOOLS` + dispatch (thirteen tools now; updated
`agent/prompts.py` accordingly, including a principle-2-adjacent note in tool #3's description
pointing at exactly what kind of claim this produces). Added a 5th checklist category,
`"regulon_target"`, to `checklist_generator.py` (using a `"TF->TARGET"` convention for the
single-string `gene_or_motif` field, e.g. `"SPI1->CD14"`, so the existing per-`(category, gene)`
dedup logic needs no change) and `report_generator.py`'s label/category maps. Added one sentence
to `novelty.py`'s prompt naming `regulon_inference` alongside `peak_to_gene_links` as a way to
avoid every novel-finding candidate collapsing into "TF vs. its own motif."

**Tests (`tests/test_regulon_inference.py`):** `_overlap_pairs`/`_genes_near_peaks` on synthetic
coordinates (4 tests), plus a fully synthetic end-to-end `infer_regulon_targets` test built around
4 candidate genes designed to isolate exactly one variable each: TARGETA qualifies via promoter
motif evidence only; TARGETB qualifies via a genuinely distal peak whose accessibility really
does track TARGETB's expression; **TARGETC has a TF-motif peak just as close, but that peak's
accessibility does NOT track TARGETC's expression -- and TARGETC's own TF1-RNA correlation is
made IDENTICAL in strength to TARGETA/TARGETB's** specifically so the test proves exclusion comes
from the motif/chromatin evidence gate, not from weaker RNA correlation; TARGETD sits near no
TF-motif peak at all. `load_protein_coding_gene_coords` monkeypatched to a small fabricated table
(no real GENCODE file needed for this logic test). All 9 tests passed on the first real run against
this fixture -- the hand-computed expected distances/overlaps (worked out on paper before writing
assertions, e.g. peak_Q-to-TARGETB = 9900bp) matched exactly. 4 more tests for missing-TF/
missing-motif-data/no-match/menu-wrapper error paths.

## Step 9: Real end-to-end verification against the real (rebuilt) PBMC cache -- and a real
compute-cost mistake caught and fixed along the way

First mistake, caught by the owner: to get `motif_match` into the agent's cached fixed-core
h5mu, I initially called `get_agent_fixed_core(force_reload=True)`, which reran the ENTIRE
pipeline from scratch -- including chromVAR's slow permutation step, which the owner correctly
pointed out `_persist_motif_match` never touches (it only needs `pc.match_motif`'s output,
computed BEFORE `pc.get_bg_peaks`/`pc.compute_deviations` even run). Fixed by writing a targeted
script that loads the EXISTING cache (already-correct `chromvar_deviations`, gene activity, QC,
clustering, all untouched) and reruns ONLY peak-sequence lookup -> GC bias -> JASPAR motif
matching -> persist, then writes back. Real timing: ~2 minutes instead of the ~30-60+ minutes a
full rebuild would have cost. Same fix applied to the shareseq cache (133,743 peaks x 879 motifs
persisted, ~2 minutes).

**Second real problem, found only by actually running it on real data:** `regulon_inference(mdata,
"SPI1")` hung for 10+ minutes without finishing. Root cause, found by `cProfile` on a real run
rather than guessed: SPI1's JASPAR motif ("MA0080.7.Spi1", an ETS-family motif) matches 48,557 of
107,385 real PBMC peaks -- 45% of the whole genome's worth of accessible chromatin, a real
property of this motif's short/degenerate consensus, not a bug. The profile showed 499 of 624
seconds inside `scipy.sparse._sparsetools.get_csr_submatrix`: `rna[:, gene].X` (AnnData's own
column indexing) was being called once per unique candidate gene (12,432+ times), and slicing a
single column out of a CSR matrix is O(total nnz) per call, with AnnData's view/
`_remove_unused_categories` machinery adding further overhead on top. Fixed at the real source:
convert `rna.X`/`atac_counts` to CSC ONCE per call and index the raw matrix directly, bypassing
AnnData's `__getitem__` entirely -- cut SPI1's real runtime from 624s to 43s (verified: identical
result, 12,432 candidates / 1,134 significant, both before and after). Also fixed, while profiling
exposed it: the per-pair Spearman correlation loop was replaced with a vectorized rank-transform +
chunked dot-product (`_rank_zscore_rows`/`_spearman_pvalue_from_rho`, numerically verified
identical to `scipy.stats.spearmanr` on real test data before trusting it), and a transparent,
logged cap (`DEFAULT_MAX_DISTAL_TF_MOTIF_PEAKS=2000`) was added for the DISTAL search specifically
when a motif is this abundant -- promoter evidence stays exhaustive regardless, since it's cheap
(vectorized overlap only, no per-pair correlation).

**Real results, PBMC** (SPI1/GATA3/PAX5) and **shareseq** (FOXA1/FLI1 on the real T-47D/A-673
cell lines this dataset actually contains, resolved via `detect_identity_columns` -- free, no LLM;
NKX2-1 skipped, not present in this dataset's gene symbols): see the owner's literature-check
conversation turn for the full table. Two independently-read PubMed confirmations (not title-only):
SPI1->HCK (PMID 21993313, "PU.1 was capable of transactivating the HCK promoter... ChIP...
provided evidence that PU.1... bound to the HCK promoter in vivo") and GATA3->BCL11B (PMID
35705252, "GATA-3 induces Bcl11b"). FLI1 found 0 significant targets on A-673 -- plausibly because
A-673 is driven by the EWSR1-FLI1 fusion, which binds GGAA-microsatellite repeats differently than
wild-type FLI1's canonical JASPAR motif, a real reportable limitation rather than a tool failure.

**No API calls were made anywhere in this step** -- the profiling, the perf fix, and the real runs
above are all local compute; the literature check used only free NCBI E-utilities calls.

## Step 10: Cost-saving shortlist-first path for checklist generation + novelty judging

Following a discussion with the owner about why `generate_checklist`'s `run_agent` calls cost
real money when the PubMed calls inside them are free: the cost was never the PubMed step -- it's
the open-ended, multi-turn EXPLORATION (discovering dataset identity from scratch, brainstorming
candidates with no hint, trying up to 4 per category, reading abstracts and judging fit). When a
dataset's real identity is ALREADY retrievable for free from its own metadata
(`detect_identity_columns`, no LLM -- true for shareseq's real cell-line names, false for the
public PBMC dataset, which carries none), most of that exploration is redundant.

New module `agent/shortlist.py`:
- `resolve_cheap_identity(mdata)` -- wraps `detect_identity_columns`, prefers a human-readable
  `cell_line_name` column over ID-coded ones (ACH-XXXXXX), defensive (returns `None` rather than
  raising on malformed/bare `mdata`, so this cost optimization can never be the reason the
  pipeline crashes).
- `generate_shortlist_candidates(identity, category)` -- ONE-SHOT, tool-free recall on
  `SHORTLIST_MODEL` (Haiku, the cheapest tier in `agent.loop.MODEL_COSTS`) -- no exploration, just
  "name up to 5 well-known candidates," forced-tool-choice for structured output.
- `verify_candidate_in_data(mdata, category, primary)` -- ZERO-LLM dispatch to the same
  already-free deterministic analysis-menu functions (`tf_motif_correlation`, `peak_to_gene_links`,
  `regulon_inference`, `marker_is_recovered`, `best_motif_match` against `chromvar_motif_names`)
  the expensive path would call anyway, just invoked directly instead of through a model's
  tool-use turn.
- `shortlist_checklist_items_for_category` -- orchestrates: recall -> verify in data (free) -> ONE
  targeted `search_pubmed`+`fetch_pubmed_abstracts` (free) -> ONE cheap single-abstract judgment
  (`judge_abstract_supports_claim`, Haiku) -> up to `needed` confirmed `ChecklistItem`s. Returns
  fewer than `needed` (even zero) honestly if the shortlist doesn't pan out -- never forces a weak
  candidate.
- `cheap_already_known_check(finding)` -- for novelty: ONE free search + one cheap judgment;
  if a proposed "novel" finding is actually already well-established, strikes it down right there
  -- the expensive adversarial Judger agent call is skipped entirely for that finding, since
  "already known" is already settled.

**`checklist_generator.generate_checklist`** gained `categories`/`existing_items`/`use_shortlist`/
`shortlist_model` parameters. Per category: try the shortlist path (needed=3); only categories
that don't reach 3 get a SCOPED fallback `run_agent` call (one call covering just the still-short
categories, with the already-resolved identity passed in as a hint so the agent doesn't re-derive
it) -- `categories` not passed at all default to the full original behavior unchanged.
`existing_items` for categories NOT in `categories` pass through completely untouched, verbatim,
never re-verified or re-cited -- the mechanism the owner asked for to avoid re-running/re-paying
for RNA-marker/motif categories when only adding the two new regulon-era categories to an
already-generated report. If no identity is cheaply available at all (the public PBMC case),
`resolve_cheap_identity` returns `None` and behavior is byte-for-byte identical to before this
existed -- zero risk of regression for that dataset.

**`agent/novelty.py`**: `_novelty_question` now states an explicit PRIORITY ORDER -- the first 3
candidate attempts must be `regulon_inference`- or `peak_to_gene_links`-based (TF-target-gene or
peak-to-gene claims), never `tf_motif_correlation`-based (TF-vs-own-motif); only after those first
3 attempts fail to survive their own literature self-check may a TF-vs-own-motif candidate be
tried for the remaining attempts. `judge_novel_findings` gained `use_cheap_prefilter`/
`shortlist_model`, trying `cheap_already_known_check` before ever spending a real adversarial
Judger call; `check_novelty_negative_control` passes both through.

**Deliberately NOT implemented** (explicit owner instruction): provenance tagging of which path
(cheap shortlist vs. expensive exploration) produced a given item -- the report stays silent on
that distinction. Also not implemented: fine-grained token/cost tracking for the shortlist's own
(small) Haiku calls merged into the returned `AgentRunResult` -- noted here as a real, known gap,
not silently dropped; the zero-cost placeholder result returned when the shortlist alone suffices
under-reports the true (still small) cost of the Haiku calls that ran.

**Tests** (`tests/test_shortlist.py`, 17 tests; plus 6 new in `test_checklist_generator.py`, 3 new
in `test_novelty.py`): every LLM call (`_call_shortlist_model`) and every network call
(`search_pubmed`/`fetch_pubmed_abstracts`) is mocked -- confirmed via a 30s-timeout test run that
nothing in this module makes a real network or API call during tests. All existing tests in both
files pass UNCHANGED (0 modifications) for `checklist_generator.py`'s suite -- `resolve_cheap_identity`'s
defensiveness on `mdata=object()` means the old tests silently skip the new shortlist branch
exactly as before. `test_novelty.py`'s 4 pre-existing `judge_novel_findings`/
`check_novelty_negative_control` calls DID need one added keyword (`use_cheap_prefilter=False`)
each, since there's no equivalent "automatically skip" signal for the novelty prefilter (it
depends only on finding text, not on `mdata`) -- a real, acknowledged test-only change, not a
behavior change.

**No API calls were made building or testing this** -- the actual cost-saving effect (real Haiku
shortlist calls against the real PBMC/shareseq data, replacing real Sonnet/Opus exploration) is
applied separately, with the owner's explicit go-ahead, directly against the two committed example
reports -- see the owner's conversation turn authorizing it and whatever follow-up entry is added
once that run completes.

## Step 11: Retired "tf_motif_tracking" as a checklist category entirely

Owner's follow-up, directly motivated by everything above: "TF's own RNA vs. its own motif"
shouldn't be a known-biology checklist category any more at all -- it's exactly the narrow,
easily-already-known pattern this whole phase has been working to move past, and now that
"regulon_target" (TF vs. a SPECIFIC OTHER gene it regulates, real motif+chromatin evidence) exists
as a structurally stronger replacement, there's no reason to keep the weaker one as a checklist
category too. `CHECKLIST_CATEGORIES` is now `("rna_marker", "motif", "peak_to_gene",
"regulon_target")` -- four categories, not five. `tf_motif_correlation` (the underlying tool)
is untouched and still real/usable elsewhere (e.g. as one of several verification tools, or inside
novelty's TF-vs-own-motif fallback after 3 failed regulon/peak-to-gene attempts) -- only its
checklist-category status was retired.

Updated: `CHECKLIST_QUESTION`'s numbered list (four categories now, with an explicit note on why
`tf_motif_tracking` was dropped) and verification-tool mention; `_FALLBACK_CATEGORY_DESCRIPTIONS`
and `agent/shortlist.py`'s `_CATEGORY_INSTRUCTIONS`/`verify_candidate_in_data` (removed the now-
dead `tf_motif_tracking` branch and its now-unused `tf_motif_correlation` import);
`report_generator.py`'s `_render_section3` labels/`by_category` dicts. Existing tests that used
"tf_motif_tracking" purely as an arbitrary example category string (dedup/rendering logic tests,
which don't validate against `CHECKLIST_CATEGORIES` at all) were updated to use a real current
category (`peak_to_gene`/`regulon_target`) instead, to avoid the test suite itself implying the
retired category still exists; the one test that asserted the exact category tuple was updated to
the new 4-tuple. One dead test (`verify_candidate_in_data` dispatching on `tf_motif_tracking`) was
deleted outright rather than kept as a test for unreachable code.

All 60 tests across `test_checklist_generator.py`/`test_novelty.py`/`test_shortlist.py`/
`test_report_generator.py` pass; full-suite collection (242 tests) is clean. Confirmed no other
`src/` file referenced `tf_motif_tracking` or `CHECKLIST_CATEGORIES` outside `checklist_generator.py`
itself (`grep` across `src/`), so no other code needed updating.

## Step 12: Updating the two committed example reports with real API calls

With the owner's explicit go-ahead: adding "Distal peak-to-gene links" and "TF-target-gene
regulon links" subsections to section 3 of both `reports/examples/tenx-cell-ranger-report.md`
(PBMC) and `reports/examples/shareseq-multi-cell-lines-report.md`, without re-running or
re-rendering sections 1/2/4 or section 3's RNA-marker/motif subsections (preserved verbatim),
regenerating section 5 (novel findings) with the new TF-target-priority ordering, and leaving
each report's cost header/fault-injection cost table untouched per explicit instruction.
Moderate-price model (`claude-sonnet-5`) for exploration/novelty; the shortlist step itself still
uses Haiku. Negative-control subsection left untouched -- shuffled data has no real signal to
find regardless of which analysis type is tried first, so re-running it wouldn't demonstrate
anything different, only cost money.

**Real paid run #1 failed -- a genuine production bug, not a test artifact.** The PBMC run (no
cheap identity available, so straight to full `run_agent` exploration) called `regulon_inference`
on EOMES and RUNX3 -- both abundant motifs (16,584/17,814 of 107,385 peaks respectively) -- and
each returned THOUSANDS of candidate targets (7,654 / 7,771) with NO cap on the returned `targets`
list at all (unlike `peak_to_gene_links`, which already had a `top_n` cap -- an oversight in the
original `infer_regulon_targets` implementation that Step 9's abundant-motif fix never addressed,
because that fix was about COMPUTE time, not OUTPUT size). Both full lists got serialized into
their tool-result messages and sent back to the model; by the next turn the conversation's total
prompt hit 2,281,949 tokens and the API rejected it outright (`400: prompt is too long... >
1000000 maximum`), crashing the run after ~15 real turns of real spend with zero usable output.

**Fixed properly, not just patched around:** added `top_n` (default 30) and `target_gene`
parameters to `infer_regulon_targets`. The returned `targets` list is now ALL significant targets
first (up to `top_n`), then the strongest non-significant ones to fill out context -- but
`n_candidate_targets`/`n_significant_targets`/`any_significant_target` still reflect the TRUE full
counts, and the underlying correlation + BH-correction still runs over every real candidate
regardless of `top_n` (only the returned LIST is bounded, not the statistics). Critically, capping
the list creates a real risk for the checklist/novelty use case specifically: verifying one named
candidate (e.g. "does SPI1 regulate CD14") could get a false negative if that gene is significant
but just outside the top-N cut. Fixed by adding `target_gene`: when passed, that gene's real entry
is guaranteed to appear in `targets` even if it wouldn't otherwise make the cut. Threaded through
`menu/regulon_inference.py`, the agent-facing tool schema + dispatch in `agent/loop.py` (so the
model itself can pass `target_gene` once it has a specific candidate in mind), and
`agent/shortlist.py`'s `verify_candidate_in_data` (which now always passes `target_gene=target`
for the `regulon_target` category). Updated `CHECKLIST_QUESTION`/`_FALLBACK_CATEGORY_DESCRIPTIONS`/
`CORE_SYSTEM_PROMPT`'s tool #3 description to explicitly instruct the agent to always pass
`target_gene` once it has a specific candidate, not just `tf_gene` alone.

Verified the fix directly against the exact real cases that caused the failure (free, local, no
API cost): EOMES and RUNX3 against the real PBMC cache now return `len(targets)==30` (down from
7,654/7,771) at ~10KB serialized (down from whatever blew past 1,000,000 tokens) -- true counts
(`n_candidate_targets`, `n_significant_targets`) unchanged. Added
`test_infer_regulon_targets_truncates_but_target_gene_always_included` (a 6-target synthetic
fixture; the expected strength ranking is read from an UNCAPPED call's own output rather than
hand-derived from the cyclic-shift construction, since that correlation isn't simply monotonic in
shift size) confirming: the true count is cap-independent, the weakest candidate is excluded under
a small `top_n`, and reappears when requested via `target_gene`. All 84 tests across the directly
affected files pass.

**Real paid run #2, tenx (PBMC) -- completed.** `generate_checklist(categories=("peak_to_gene",
"regulon_target"))` ($8.82, 40 turns) found 2 confirmed peak_to_gene items (CD8A, CD69) + 1
rejected (IL2RA) and 3 confirmed regulon_target items (PAX5->CD19, EOMES->PRF1, KLF4->CD14);
`propose_novel_findings`/`judge_novel_findings` ($0.83 propose + $0.42+$0.18+$0.29 judging) found
3 findings (CEBPB->FGL2 struck down as a monocyte-lineage confound; MS4A1 distal peak struck down
as a cell-type-identity confound spanning the whole locus; TCF7->IL7R struck down as already-known
co-expression with the opposite causal direction in one retrieved paper) -- real tenx/CEBPB/MS4A1/
TCF7 numbers all cited directly in `reports/examples/tenx-cell-ranger-report.md`'s updated section
3/5. `tf_motif_tracking` retired as a checklist category entirely per Step 11 (the old SPI1/MAFB/
CEBPA confirmed items and JUNB/KLF4 rejections were removed from the report, not just appended to).

**Real paid run #2, shareseq -- completed in a separate pass**, after the owner's follow-up
clarified (with a filesystem-full/lost-context episode in between, recovered via this doc + the
session's own tool-call history) that the shareseq half of this update had been paused, not
abandoned, because the committed `shareseq-multi-cell-lines-report.md` had accidentally been
overwritten with tenx content at some point before this session, and that was separately fixed
by the owner first. Same fix (`shortlist.py`'s `peak_to_gene`/`regulon_target` candidate-format
prompts) applied before retrying: an earlier attempt at this exact shareseq checklist call had
Haiku return a malformed `primary` string for `peak_to_gene` (a compound
"CELL_LINE_DISTAL_ENHANCER->MYC" instead of a bare gene symbol) and a bare TF name with no "->"
for `regulon_target`, both of which fail `verify_candidate_in_data` silently and force an
unnecessary (small, ~3-turn, killed-on-timeout) fallback -- fixed by making `_CATEGORY_INSTRUCTIONS`
state the required format explicitly with real examples, plus a defensive `rsplit("->", 1)`
normalization for `peak_to_gene` specifically (no equivalent safe recovery exists for a
`regulon_target` candidate missing its target entirely). Retried clean: checklist ($2.55, 20
turns) found 0 confirmed items in either new category but 6 real, literature-grounded REJECTIONS
(ELF3/EHF/TGIF1 peak_to_gene; ELF3->EHF/NFE2L2->NQO1/SPDEF->MUC5AC regulon_target -- an LUAD core
regulatory circuit and two other well-established relationships that do NOT replicate in this
specific pooled 8-cell-line dataset's own data); novelty ($0.38 propose + $1.22+$1.16+$2.54
judging) found MYCN->GNG7 (struck down as a cross-cell-line cluster confound), FOXA1->KYNU (the
Judger call ended without calling `record_judger_verdict` even after the automatic nudge -- a
known, previously-documented failure mode, reported as "NO VERDICT RECORDED" rather than silently
dropped or re-run at further cost), and EGFR-no-distal-link (a `no_signal_or_concern` finding that
still went through the full Judger under this pipeline's current behavior -- struck down as a
pooled-dataset scoping confound). `tf_motif_tracking`'s old FOXA1 confirmed item and STAT3/GATA3/
FOS rejections were likewise removed from the shareseq report, matching tenx's treatment exactly.

**Follow-up correction, same session: real numbers + a lowered threshold for shareseq.** The
owner asked for the actual rho/p/q numbers behind every peak_to_gene/regulon_target confirmed AND
rejected item in BOTH reports (not just prose claims), pulled from the already-paid-for logs
rather than a re-run. Recovered from `logs/agent_runs/20261006T004747Z.md` (tenx) and
`logs/agent_runs/20261006T040207Z.md` (shareseq, which conveniently contains the model's own
summary table with most of the real numbers already in it) directly; the handful the model only
reported as a rounded "≈" figure, or that the 1000-char log truncation cut off entirely (PAX5->
CD19, KLF4->CD14's exact values; EHF's exact rho, right at the position a 0.1-vs-0.2 threshold
decision would hinge on) were recovered via a genuinely free fallback -- calling
`peak_to_gene_links`/`regulon_inference` directly against the same already-cached, unchanged data
(deterministic, $0 API cost, not an agent run) rather than leaving them as approximations.
Separately, per the owner's instruction, shareseq's (not tenx's) `peak_to_gene`/`regulon_target`
confirmation bar was lowered to |rho|>=0.1 (q<0.05 still required) -- applying this to the
recovered exact numbers flips **NFE2L2->NQO1 from rejected to confirmed** (rho=0.1319, q=5.2e-22);
ELF3/EHF/TGIF1 (rho=0.052-0.098) and ELF3->EHF (rho=0.091) all stay rejected even at the lowered
bar; SPDEF->MUC5AC stays rejected for a more fundamental reason (never even a motif-supported
candidate, independent of any threshold). Both reports' "Confirmed by data"/"Rejected by data"
subsections for these two categories now carry the real numbers inline, and shareseq's carries an
explicit note on the lowered, dataset-specific bar so a reader isn't confused by a number that
wouldn't clear the tool's own default 0.2 cutoff. The "Cost breakdown" section at the bottom of
both reports was left completely untouched (matches the ALREADY-established precedent from the
first tenx update -- confirmed via `git diff` that section genuinely wasn't touched then either),
since it documents the original full `generate_report` run's cost, not a running total of every
later patch.

## Step 13: Verified `make demo` itself still works end-to-end, unrelated to any of the above

Before touching the demo cache, the owner asked to actually verify (not assume) that the
README's documented quickstart -- `conda env create`, `cp .env.example .env`, `make demo` --
genuinely still works after everything above, specifically on questions that DON'T touch
`peak_to_gene_links`/`regulon_inference` (since those were already known, at this point, to be
broken in the demo specifically -- see Step 14). Ran the literal default demo question plus two
more (`"does CD3E mark T cells in this dataset?"`, the exact example from the script's own
docstring; `"is GATA3 a marker of any cell cluster in this dataset?"`, a genuine negative result
the agent correctly reported as an absence rather than forcing). All three: real runs, no crashes,
no tool errors, 3-4 turns, $0.009-$0.014 each, 7-10s each -- comfortably inside the README's
"~10-25s, ~$0.02" claim. Confirms none of this phase's prompt/tool-list/checklist-category changes
regressed the demo's core promise for the large majority of questions that don't touch the two
newest tools.

## Step 14: Fixed the demo cache itself so `peak_to_gene_links`/`regulon_inference` work there too

Real gap, found while answering the owner's direct question ("does the demo also have regulon
analysis") rather than assumed: `agent/demo_fixed_core_cache.py`'s `_trim_to_tool_readable_fields`
predates these two tools entirely and stripped every one of the four fields they need --
`atac.var` wiped to zero columns (no `chrom`/`start`/`end` to locate peaks genomically),
`atac.varm` cleared entirely (no `motif_match`), `atac.uns`'s keep-list excluded
`motif_match_names`, and `atac.layers["counts"]` dropped with `.X` zeroed (no real accessibility
to correlate against). Both tools would fail -- gracefully (a caught exception/explicit error
dict, not a crash), but fail -- on any demo question that tried to use them. The default demo
question (SPI1-vs-own-motif) never exercises either tool, so this went unnoticed until directly
asked about.

**Fixed the source function** (so a FUTURE from-scratch `build_demo_cache()` rebuild stays
correct) to keep `chrom`/`start`/`end`, the real `counts` layer, and `motif_match`/
`motif_match_names` -- the last two filtered down to only motifs whose real TF-name token
(`core.motif_deviations.motif_name_tokens`, factored out of `_motif_match_rank` for this reuse)
matches a gene actually present in the demo's own RNA panel. This filter is lossless for real
functionality: `regulon_inference` already rejects any OTHER TF at its very first gate (`tf_gene
not in rna.var_names`) regardless of whether its motif column exists, so those columns were
already 100% dead weight for this specific RNA panel -- confirmed 468 of JASPAR's 879 motifs fall
in this category (non-human factors, composite-only entries, genes filtered out of this panel),
leaving 411 genuinely reachable ones for the committed demo's RNA panel specifically.

**Repaired the already-committed artifact without paying the ~17-20 min from-scratch rebuild
cost again**: wrote `scripts/backfill_demo_regulon_fields.py`, a one-time migration that reindexes
the four missing fields directly from the already-built, already-`motif_match`-persisted real
agent fixed-core cache (`data/processed/agent_fixed_core_n11909_seed0.h5mu`) -- verified first
(not assumed) that this is safe: the demo's 150 cell barcodes are a real subset of the full
cache's 11,909 (checked via actual set intersection, since a shared RNG seed at different draw
sizes gives no guarantee of a nested draw), and the demo's own independently-re-filtered-for-150-
cells 61,474-peak set is a real name-subset of the full cache's 107,385 peaks (also checked, not
assumed).

**Real, measured file-size problem caught before committing, not after**: the first version of
this backfill (before the motif-reachability filter) grew the committed demo cache from 48.5MB to
114.8MB -- over this project's own documented 100MB GitHub hard-limit (the exact ceiling the
module's docstring warns about for the FULL, untrimmed object). Diagnosed by measuring the actual
h5 dataset sizes directly (not guessed): `varm/motif_match`'s sparse CSR storage alone cost 58.3MB
-- MORE than storing it dense would, because real density here is ~21.5% and sparse CSR's int32
index overhead stops paying for itself well before that density. The reachability filter above
(879 -> 411 motif columns) fixed this as a side effect of being the correct thing to do anyway --
final committed size: 84.8MB, comfortably under the limit, with zero loss of real functionality
(verified: every dropped motif was already unreachable for this RNA panel regardless).

**Tests**: added `test_demo_cache_peak_to_gene_links_still_works` / `test_demo_cache_
regulon_inference_still_works` (real calls against the repaired committed cache) and extended
`test_demo_cache_has_every_field_the_agent_tools_read` with assertions on the four fields. All 7
tests in `tests/test_demo_fixed_core_cache.py` pass; full-suite collection (245 tests) is clean.

**Verified against the real, repaired, committed demo cache via `make demo` itself** (not just
unit tests): `make demo QUESTION="does SPI1 regulate CD14 in this dataset, beyond just its own
motif?"` -- real `regulon_inference` call, found CD14 as a significant SPI1 target (rho=0.42,
both promoter AND a real distal-peak evidence entry), 2 turns, $0.017, ~15s.
`make demo QUESTION="is there a distal peak near CD14 whose accessibility tracks CD14
expression?"` -- real `peak_to_gene_links` call, found 6 significant distal links (best rho=0.34),
2 turns, $0.00535, ~9s. Both clean, no errors, numbers consistent with what the same tools report
against the full-scale data.

Updated `README.md`'s demo-cache size claim (46MB -> ~85MB, the real current committed size) and
added one sentence noting the demo now genuinely exercises all 13 tools, with a worked
`peak_to_gene_links`/`regulon_inference` example question, not just the default TF-motif one.

**Spot-checked the newly-fixed demo cache against 5 more real regulon-related questions** (not
just the two from the verification above) before changing anything further: `regulon_inference`
with no `target_gene` given (open-ended -- "does GATA3 regulate any specific target gene?" ->
correctly surfaced 20 real significant targets, PRR5/PBXIP1/LRP10/... with real rho/q); the exact
same question again as a sanity check on reproducibility; a genuinely negative case ("does TCF7
regulate IL7R?" -> IL7R correctly absent from TCF7's significant target list at this 150-cell
scale, and the agent correctly pivoted to a real literature search rather than forcing a claim,
surfacing PMID 15365098's actual OPPOSITE-direction mechanism -- a good demonstration of "report
real rejections honestly" working end-to-end in the demo, not just the full pipeline); and a
second `peak_to_gene_links` case (CD3E's own distal enhancer, found real, rho=0.31, 82kb away).
All 5 real, cheap ($0.005-$0.032 each), clean, no errors.

## Step 15: Changed the demo's default question to showcase `regulon_inference`

Per the owner's request: `scripts/demo.py`'s `DEFAULT_QUESTION` changed from the original
SPI1-vs-own-motif question (which only ever exercised `tf_motif_correlation`, predating
Phase 2 entirely) to the GATA3 open-target regulon question from the spot-check above --
so `make demo` with NO arguments now showcases the newer, structurally-different
`regulon_inference` capability by default, not just the narrower TF-vs-own-motif check. Verified
the new default actually runs clean (`make demo`, no QUESTION override): 2 turns, $0.00852,
~14s, same 20-significant-target GATA3 result as the spot-check. Updated the script's own
docstring and `README.md`'s demo-cache paragraph (removed the stale "24s/$0.022" real-run number,
which was from the OLD default's own call, and the stale "example" question it used) to match.
