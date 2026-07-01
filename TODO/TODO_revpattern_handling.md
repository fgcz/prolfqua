# Decoy and contaminant handling in prolfqua core (rev/contaminant patterns)

**Date:** 2026-06-30 (design); **updated 2026-07-01** (core implemented + single-path plan settled).
**Status:** **prolfqua-core machinery IMPLEMENTED**; the prolfquapp single-explicit-filtering-path
refactor is **designed + source-verified, ready to implement** (see "Single explicit filtering
path" below). No longer deferred — a synthetic decoy+contaminant fixture removes the earlier
"can't validate without decoy data" blocker.
**Driver:** the prolfquapp redesign in
`prolfquapp/TODO/TODO_protect_prolfquapp_against_target_decoy.md` (handling target+decoy FASTAs,
guaranteeing unique protein IDs). This note records the `prolfqua`-side follow-up so it is not
lost.

**Key decision update (2026-07-01):** **contaminants are KEPT and LABELLED, not removed.** This
supersedes the earlier §4 "default remove". Decoys and contaminants are now asymmetric in
*visibility*, not in whether they pass through:

- **decoys** — kept through normalization → dropped at the *fit* → NA on export → **invisible**
  (machine artifacts, must never be findings);
- **contaminants** — kept **everywhere** (normalization, fit, contrasts, export), carrying an
  `is_contaminant` / annotation `CON` flag so figures can **label** them (real proteins, shown but
  marked). No contaminant removal step; `LFQData$remove_contaminants()` is therefore dropped.

## Implementation status (2026-07-01)

**Landed in `prolfqua` (committed on `main`):**

- **Shared detectors** (`706e08af`, `R/utilities.R`): `is_decoy()`, `is_contaminant()`,
  `effective_decoy_pattern()`, `effective_contaminant_pattern()` — one database-agnostic
  implementation (config pattern ∪ built-in defaults; empty/`NULL`/`"a^"` → defaults only).
  Resolves **F4** (effective-pattern exposure). Tests: `test-decoy-contaminant-detect.R`.
- **Config slots** (`4166581f`, `R/AnalysisConfiguration.R`): optional `pattern_decoys` /
  `pattern_contaminants` (default `NULL` = off; round-trips through `R6_extract_values`).
- **LFQData methods** (`4166581f`, `R/LFQData.R`): `remove_decoys()` (targets-only fit),
  `decoy_proportion()` / `contaminant_proportion()` (QC, per modelling-level `subject_id`; decoy
  status derived from the top-level protein id — **F5**). Tests: `test-decoy-contaminant-lfqdata.R`.
  _(An earlier `remove_contaminants()` was removed in `b85b535e` — contaminants are kept+labelled,
  see the top-of-file decision.)_
- **Targets-only fit** (`4166581f`, `R/build_contrast_analysis.R`): when `pattern_decoys` is set,
  decoys are dropped before the fit, so the fit + the shared variance pool (limma prior / DEqMS
  variance-count trend) see targets only — **F1**. Opt-in: `NULL` leaves existing behaviour
  untouched. Tests: `test-build-contrast-decoy-drop.R`. Full suite: 1047 pass, 0 fail.

**Landed in `prolfquapp` (committed on `master`):**

- **Detector consolidation** (`e3a2bcd`, `R/R6_ProteinAnnotation.R`): `.detect_decoy_ids` is now a
  thin wrapper over `prolfqua::is_decoy` — the byte-for-byte duplicate default-prefix set is gone,
  so annotation de-duplication and the prolfqua quant layer share one detector. Behaviour
  identical; ProteinAnnotation tests: 26 pass. WU347806 end-to-end re-verified (SE builds, 4029
  unique proteins).

## Single explicit filtering path — investigation + plan (2026-07-01)

Source-verified by a fan-out investigation across all readers + the export path (workflow
`wf_9b95c4fe-aa1`, 6 agents). The goal (your directive): **one** explicit quant filtering path — no
alternative/redundant filtering. The annotation↔quant join must **preserve** all quant rows, never
filter them.

### Load-bearing findings (verified against source)

- **GAP A — decoy-drop-at-fit is dead code in prolfquapp.** The only decoy drop before the fit is
  the core guard `build_contrast_analysis.R:105-119` (`if (!is.null(cfg$pattern_decoys))
  lfqdata <- lfqdata$remove_decoys()`). prolfquapp **never calls** `build_contrast_analysis()`;
  `DEAnalyse$build_facade` (`R6_DEAnalyse.R:168`) constructs `facade_class$new(self$lfq_data, ...)`
  directly, and no facade `$new()` does decoy handling. So the fit-drop guard is unreachable from
  prolfquapp.
- **GAP B — the LFQData config never carries the patterns.** `AnalysisConfiguration` defaults
  `pattern_decoys`/`pattern_contaminants` to `NULL`; **no** prolfquapp preprocessor assigns onto
  `lfqdata$get_config()$pattern_decoys` — patterns are passed only to `ProteinAnnotation$new()`.
  So even if Gap A were reachable, the gate `!is.null(cfg$pattern_decoys)` is always false.
- **The inner-join IS the hidden/redundant filter.** `ProteinDataPrep$remove_cont_decoy()`
  (`R6_ProteinDataPrep.R:67-72`) does `get_subset(clean())`; `LFQData$get_subset()` is an
  `inner_join` on `protein_Id`. Because the annotation is decoy-free at construction, this inner
  join **drops decoy quant rows at the peptide level, before aggregation** — the exact
  "alternative/redundant filtering path" to eliminate. Redundancies enumerated: decoys dropped at
  BOTH `remove_cont_decoy` and the (dead) fit guard; contaminants filtered at `remove_cont_decoy`,
  `run_dea`, and again for IBAQ (`cmd_helpers.R:458-460`); IBAQ also inner-joins the annotation
  (`aggregation_IBAQ.R:83-87`).
- **Export already does the right thing — for free.** The SE row set comes from `lfq_data_raw`
  (`R6_DEAReportGenerator.R:573-575`), and `.join_annotation` is a **`right_join`**
  (`report_helpers.R:42-49`), so if decoys survive in `lfq_data_raw` but are dropped from the fit,
  the SE/results get decoy rows with **NA** contrast stats (base-`[` NA-padding at
  `R6_DEAReportGenerator.R:624`) — automatically absent from ORA/GSEA/.rnk/volcano. **No export
  code changes needed** (fix for **F2**).
- **Quant-side prefix retention per reader** (decides whether `is_decoy`/`is_contaminant` can act on
  the quant `protein_Id`): all readers **retain both `REV_` and `CON/zz`** on the quant id
  (MaxQuant `leading.razor.protein`, FragPipe `Protein`, BGS `PG.ProteinGroups`, MSstats
  `ProteinName`, PTM readers map to `fasta.id`) **except DIA-NN**, which strips the `zz|` wrapper
  (`preprocess_DIANN.R:91`) but keeps `REV_`. → **decoys detectable on all readers' quant;
  contaminants NOT detectable on DIA-NN quant.**
- **BLOCKER C (contaminant detection on DIA-NN quant) — sidestepped by keep+label.** Since
  contaminants are no longer *removed* on the quant side, we never need to detect them there. For
  **labelling**, drive it off the annotation's `CON` flag carried onto results via the existing
  `right_join` (the annotation's `full_id`=`fasta.id` retains `zz`, so it works for DIA-NN too).
  This avoids the risky "stop stripping the DIA-NN prefix" change to the join key entirely.

### The ONE filtering path (target)

```
reader → LFQData(peptide) + ProteinAnnotation(decoy-free; CON-flagged, NOT contaminant-removed)
   │  [FIX B] stamp config$pattern_decoys (+ pattern_contaminants) onto the LFQData config
   ▼
ProteinDataPrep$remove_cont_decoy()  → NO quant filtering (retire get_subset(clean()));
   │                                    compute decoy/contaminant proportion QC only
   ▼
aggregate() → transform_data()       (decoys AND contaminants flow through; normalized on full data)
   ▼
DEAnalyse$build_facade()
   │  [FIX A] model_lfq <- if (!is.null(cfg$pattern_decoys)) lfq_data$remove_decoys() else lfq_data
   │           facade_class$new(model_lfq, ...)     ← decoys dropped ONLY here (targets-only fit, F1)
   │           (contaminants stay in the fit — real proteins — and get real contrasts)
   ▼
export / SE  ← lfq_data_raw keeps decoys+contaminants; right_join → decoys NA stats (invisible),
                contaminants real stats + CON flag (labelled)
```

- **Decoys** — dropped exactly once, at the fit (Gap A fix), gated on `pattern_decoys` (Gap B fix).
- **Contaminants** — never removed; flagged for labelling. `remove_cont` becomes vestigial for the
  default path (keep it as a future opt-in only if a "biologist wants them gone" mode is revived).
- **Alignment** — `remove_cont_decoy`'s inner join is retired; the export `right_join` (already
  quant-preserving) is the single annotation↔quant join.

### Ordered edit list (to implement)

1. **prolfqua** — revert just `LFQData$remove_contaminants()` (keep `remove_decoys`,
   `decoy_proportion`, `contaminant_proportion`, detectors); update `test-decoy-contaminant-lfqdata.R`.
2. **prolfquapp Gap B** — stamp `pattern_decoys` (+ `pattern_contaminants`) from
   `processing_options` onto every LFQData config produced in `ProteinDataPrep` (initialize +
   after `aggregate()` + after `transform_data()`/`transform_peptide_data()`; aggregation/transform
   build fresh configs, so re-stamp — a private `.stamp_patterns(lfq)` helper). Defaults:
   `pattern_decoys="^REV"`, `pattern_contaminants="^CON|^zz"` (`R6_AppConfiguration.R:20-23`).
3. **prolfquapp** — retire `remove_cont_decoy()`'s `get_subset(clean())`: **no quant filtering**
   (contaminants kept+flagged, decoys kept); keep it only for QC logging + `decoy_proportion`.
   Ensure the annotation keeps contaminants **flagged** (`CON`) rather than dropped, and that the
   `CON` flag is carried onto results for labelling.
4. **prolfquapp Gap A** — in `DEAnalyse$build_facade`, drop decoys immediately before
   `facade_class$new(...)`, gated on `cfg$pattern_decoys`; apply to the **SAINT branch**
   (`:156-160`) and the **nested** `lfq_model` path (`cmd_helpers.R:367-393`) too. Leave
   `lfq_data`/`lfq_data_raw` untouched so export keeps decoys.
5. **prolfquapp IBAQ** — align `cmd_helpers.R:458-460` / `aggregation_IBAQ.R:83-87` to the single
   path (stop using `get_subset(clean())` as a decoy filter; IBAQ may drop decoys explicitly as an
   output-table choice, but not via the shared alignment).
6. **QC** — surface `decoy_proportion()` (empirical-FDR signal) in the summary; keep
   `contaminant_proportion()`.

### Synthetic fixture (removes the "no decoy data" blocker)

New `prolfquapp/tests/testthat/test-single-filtering-path.R`, built from
`prolfqua::sim_lfq_data_peptide_config()` + `prolfquapp::add_RevCon()` (injects ~10% `REV_`, ~5%
`zz`) + a matching decoy-free/contaminant-flagged annotation. Assertions:

- after `remove_cont_decoy`: **decoys present** in quant, contaminant handling per decision, decoy
  proportion > 0;
- after `aggregate`/`transform`: **decoys survive** normalization;
- after `build_default`/`get_annotated_contrasts`: **no decoy in contrasts** (targets-only fit),
  contaminants present in contrasts (kept), and **decoys still present in `lfq_data_raw`**;
- after `make_SummarizedExperiment`: **decoy rows present with NA contrast stats**; contaminant
  rows present with real stats + `CON` flag;
- a **standalone decoy** (no forward twin) survives to export with NA stats;
- `pattern_decoys = NULL` → decoys are NOT dropped at the fit (gate respects `NULL`).

### Adversarial risks to guard

1. **Pattern not re-stamped after aggregate/transform** → fresh cloned config has `NULL`
   `pattern_decoys` → Gap A gate fails, decoys silently re-enter the fit. Mitigation: stamp on every
   produced LFQData; log when the pattern is absent at `build_facade`.
2. **Decoys re-enter the variance pool if a branch is missed** (SAINT / nested). Mitigation: edit 4
   covers all facade-construction sites.
3. **Vendor decoy prefixes not matched by the anchored defaults** (`.default_decoy_prefixes`) →
   under-counted decoys / bad empirical FDR. Mitigation: set `processing_options$pattern_decoys` per
   reader as needed.
4. **`grepl("", x)` foot-gun** — route all detection through `is_decoy`/`is_contaminant` (guarded),
   never a raw `grepl(cfg$pattern_decoys, ...)`.
5. **Consumers assuming `nrow(SE) == nrow(contrasts)`** now break (decoy rows padded). Audited: ORA
   background reads the decoy-free `row_annot` (stays target-only, correct); the Quarto SE report
   `na.omit()`s (drops decoys, fine).

### Status — implemented 2026-07-01 (`prolfquapp` `fdbf171`, `prolfqua` `b85b535e`)

- **Done:** edit-list items 1–4 + 6. `LFQData$remove_contaminants()` deleted; Gap B stamped in
  `DEAnalyse$initialize` (+ `ProteinDataPrep$initialize` for QC); Gap A drops decoys in
  `build_facade` (targets-only fit, covers nested + SAINT); `remove_cont_decoy()`'s
  `get_subset(clean())` filter retired; `cont_decoy_summary()` reports `percentOfDecoys`.
  New `test-single-filtering-path.R` (11 assertions). prolfquapp 250 pass / prolfqua 1047 pass;
  WU347806 re-verified (SE builds, 4029 unique).
- **IBAQ (item 5) — aligned.** `cmd_helpers.R` no longer pre-filters via `get_subset(clean())`;
  IBAQ runs on a copy of the full LFQData. `compute_IBAQ_values` still inner-joins the annotation
  for `protein_length` / `nr_tryptic_peptides` (a **functional** metadata join, not a
  decoy/contaminant filter), so contaminants are kept + labelled and decoys need no explicit drop
  (per decision: "no need to drop decoys from IBAQ").
- **Postponed — figure labelling.** The `CON` flag now rides the export `right_join` onto the
  results table, but the volcano / report is **not yet** wired to visually mark contaminants
  (colour / shape / legend). This is a reporting-layer change (`ContrastsPlotter` / the Grp2
  report Rmd) deferred to a later pass. Decoys are already invisible (NA stats) and need no
  labelling.

## Background

`prolfqua`'s `LFQData` / `AnalysisConfiguration` currently has **no** decoy/rev concept at all
(the only `rev` is "reverse hierarchy order", unrelated). All decoy knowledge lives in
prolfquapp's `ProteinAnnotation`.

In the prolfquapp redesign, `ProteinAnnotation`:
- is always **decoy-free** (decoys never belong in a protein annotation),
- **detects** the decoy/rev pattern (configured pattern, else a built-in default prefix set),
- **exposes** the detected pattern via a getter (e.g. `get_rev_pattern()`).

## Design (LOCKED 2026-07-01) — quant-side rev/contaminant handling

**Status of this section:** design agreed and adversarially reviewed; **prolfqua core implemented**
(see "Implementation status"). The prolfquapp realization is planned in "Single explicit filtering
path" above. **Note:** §4 below is **superseded** by the 2026-07-01 keep+label decision (see the
top-of-file "Key decision update") — contaminants are no longer removed.

### Principle

`pattern_decoys` and `pattern_contaminants` are **not** handled through `ProteinAnnotation` for the
quant data. `ProteinAnnotation` uses the decoy pattern **only to dedup decoy entries** (already
implemented). The quant-side handling lives at the **quant layer**.

### 1. Placement: `LFQData` (option B, not preprocess)

Annotate contaminant/decoy status on the quant rows inside `LFQData` (driven by config patterns),
**not** in each preprocess function. Rationale:

- enables **display** of decoys/contaminants in volcano and QC plots (impossible if filtered early);
- **one** reusable place — all readers feed `LFQData` — instead of duplicating the logic in every
  `preprocess_*` function;
- filtering (contaminants) becomes a single clean `LFQData`/config operation.

### 2. One shared detector (promote to `prolfqua`)

Promote the id-detector (`.detect_decoy_ids`, plus a contaminant analogue) into `prolfqua` so
`ProteinAnnotation` (dedup) and `LFQData` (flagging) use **one** implementation. Detection =
**configured pattern ∪ built-in default prefixes**; empty/`NULL` → defaults only.

This resolves the crux question ("what if `ProteinAnnotation` detects a different pattern, or the
config pattern is not set?"): there is **no divergence** — both call the same detector with the same
(config ∪ defaults) semantics, so flagging is consistent and works even when the config pattern is
unset. `AnalysisConfiguration` gains optional `pattern_decoys` / `pattern_contaminants` slots
(default: none → detector falls back to defaults); prolfquapp sets them from the config.

### 3. Decoys — flag, normalize-with, drop at modelling, NA on export (revised after review)

Two levels of control:

- **prolfquapp preprocess gate — `remove_decoys` (already exists: `R6_AppConfiguration.R`).**
  Decides whether decoy protein rows enter `LFQData` at all. `TRUE` → stripped in preprocess
  (biologist default; DIA-NN reports forwards only, so usually moot). `FALSE` → decoys enter
  `LFQData`, flagged.
- **`LFQData` (when decoys are present):** `is_decoy` flag from the shared detector, then:
  - **Normalization** (vsn / robscale / quantile) runs on the **full, unfiltered** data — decoys
    included. Normalization should reflect the measured intensity distribution; the low-% shift is
    accepted (addresses review finding F3).
  - **Modelling drops decoys** — the fit + variance moderation see **targets only**. This is the
    fix for **F1**: decoys must not enter the limma prior (`squeezeVar`-style) or the DEqMS LOESS
    variance–count trend, where they would perturb *target* q-values independently of the BH
    argument.
  - **Export left-joins** the (full) abundance matrix onto the (targets-only) contrast results →
    decoy rows carry their abundances but **NA** `diff` / `p.value` / `FDR` / `statistic`. So
    decoys are present in the abundance matrix + SummarizedExperiment (abundance-based QC works)
    yet **automatically absent** from the significant set, ORA, GSEA `.rnk`, and the volcano
    (NA stats → not a finding, not plottable). This is the fix for **F2**.
- **Consequence:** since decoys are not modelled, the "decoys should be non-significant"
  FDR-validation volcano is **not** available (no decoy contrasts). Abundance-level decoy QC
  (proportion, intensity distribution) still is. If FDR-validation is wanted later, decoys need a
  *separate* fit that does not share the target variance pool.
- **QC:** report decoy proportion among quantified entries (empirical-FDR signal), from the
  abundance matrix / SE.
- **Supersedes** the interim prolfquapp behavior where `ProteinAnnotation$clean()` pattern-gated-
  removes decoys from quant — that removal moves out; decoys become preprocess-gated + LFQData-flagged.

### 4. Contaminants — keep + label (SUPERSEDED — updated 2026-07-01)

> **Superseded.** The original §4 (below, struck through in intent) had contaminants
> **pattern-gated, default remove**. The 2026-07-01 decision is **keep + label, never remove.**

**Current decision:**

- Contaminants are **real** proteins (keratin, trypsin, BSA). They are **kept everywhere**
  (normalization, fit, contrasts, export) and **flagged** so figures can label them (e.g. distinct
  colour/shape in the volcano). Keeping them in the fit is statistically fine — unlike decoys, they
  carry real variance.
- **No contaminant removal step.** `LFQData$remove_contaminants()` is dropped. Labelling is driven
  by the annotation `CON` flag (from `annotate_contaminants`, detected on the prefix-retaining
  `full_id`) carried onto results via the export `right_join` — this works for **all** readers
  including DIA-NN, sidestepping Blocker C.
- `remove_cont` (`processing_options`) is **vestigial** for the default path; retain only if a
  future opt-in "remove contaminants" (biologist) mode is revived — additive, not default.
- **QC:** still report contaminant proportion (`contaminant_proportion()`).

_Original (superseded) text: `LFQData` flagged `is_contaminant` only when `pattern_contaminants`
configured; `remove_cont` default = remove (biologist) / flip to keep+flag (bioinformatics). The
keep+label decision makes "keep+flag" the single behaviour._

### 5. What moves out of `ProteinAnnotation` when we implement this

The interim prolfquapp state routes quant contaminant + (pattern-gated) decoy removal through
`ProteinAnnotation$clean(contaminants=)` + `get_subset` (an inner join). This work **retires that
inner-join filter entirely** (`remove_cont_decoy` stops filtering quant), leaving
`ProteinAnnotation` using the pattern **only** for dedup and keeping the `CON` flag **for
labelling** (not removal). Decoy exclusion moves to the **fit** (`build_facade`); contaminants are
kept throughout. The only annotation↔quant join left is the export `right_join` (quant-preserving).

## Resolved decisions

| Question | Decision |
|---|---|
| Where — preprocess vs LFQData? | **LFQData** (enables display, centralizes) |
| Detector location? | **prolfqua** (shared by annotation + LFQData) |
| Different / unset config pattern? | one detector, config ∪ defaults → no divergence |
| Do decoys enter LFQData? | prolfquapp preprocess `remove_decoys` gate decides |
| Ever drop decoys once in LFQData? | kept + **normalized-with**; **dropped from the model fit only**; NA stats on export |
| Write decoys to SE? | **Yes** — abundance + flag (NA stats), so abundance-based SE-QC works |
| Decoys in volcano / ORA / GSEA / results? | **No** — NA stats make them absent automatically |
| Ever drop contaminants? | **No (updated 2026-07-01)** — kept everywhere + flagged for labelling; `remove_cont` vestigial |
| Contaminants in volcano / results? | **Yes** — kept in the fit, shown, **labelled** via the `CON` flag |
| Where is the decoy drop coded in prolfquapp? | `DEAnalyse$build_facade` (prolfquapp bypasses `build_contrast_analysis`) — Gap A |
| How does prolfquapp turn the machinery on? | stamp `pattern_decoys` onto the LFQData config in `ProteinDataPrep` — Gap B |
| The single annotation↔quant join? | export `right_join` (`.join_annotation`); the `remove_cont_decoy` inner join is retired |
| DIA-NN strips `zz\|` from quant id (Blocker C)? | moot — contaminant labelling uses the annotation `CON` flag, not quant-string detection |
| Decoy-proportion QC | `nr_decoys / nr distinct hierarchy keys` of the LFQData (per-peptide or per-protein by analysis level) |
| Proportion for nested (peptide+protein) analysis | **deferred** |

## History — why this was deferred, and why it no longer is

- Originally deferred: touching `prolfqua` core ripples to all dependents (prolfquapp, prophosqua,
  …), and the immediate WU347806 bug was fixed by the prolfquapp uniqueness guarantee alone.
- The core machinery was then implemented additively/opt-in (detectors, config slots, LFQData
  methods, targets-only fit) — no change to existing behaviour when `pattern_decoys` is `NULL`.
- The prolfquapp realization was briefly deferred for lack of decoy test data; a **synthetic
  decoy+contaminant fixture** (`add_RevCon` over simulated data) removes that blocker, so the
  single-explicit-filtering-path refactor is now planned + ready (see the section above).

## Adversarial review (2026-07-01) — findings & resolution

Done by hand + verified against source (the review sub-agent died on a connection error). The
revised decoy handling in §3 (normalize-with, drop-at-modelling, NA-on-export) exists specifically
to resolve F1 and F2.

| # | Sev | Finding (code-verified) | Resolution |
|---|---|---|---|
| **F1** | High | "~1% decoys don't move q-values" covers BH only; **variance moderation pools across all features** — DEqMS LOESS `loess(logvar ~ log2count)` over every protein (`prolfqua/R/ContrastsModeratedDEqMS.R:60`) and the limma prior — so decoys in the fit perturb *target* q-values. | **Drop decoys at the model fit** (targets-only), §3. |
| **F2** | Med/High | Design only de-emphasised decoys in the volcano; silent on results table / ORA / GSEA (`prolfquapp/R/R6_DEAReportGenerator.R:31` `.write_ORA`). A decoy must never be an enrichment finding. | **Left-join export → NA stats**; decoys automatically absent from significant set / ORA / GSEA / volcano, §3. |
| **F3** | Low/Med | Global normalization (vsn/quantile/robscale) is fit across the whole matrix → decoys shift it. | **Accepted deliberately** — normalize on full data; low-% shift is intended (normalization reflects the measured distribution). |
| **F4** | Med | `get_rev_pattern()` returns configured-or-`NULL`, but the shared detector adds built-in defaults → what LFQData flags (defaults) ≠ what the getter reports (`NULL`) when config is unset. | **RESOLVED:** expose the *effective* pattern (config ∪ defaults), not the raw config. |
| **F5** | Med | `is_decoy` is per-protein but `LFQData` is often peptide/precursor-level. | **RESOLVED:** flag lives at the LFQData hierarchy level (peptide rows carry `protein_Id`, so status derives from the protein). Proportion = `nr_decoys / nr distinct hierarchy keys` — per-peptide for peptide-centric, per-protein for protein-centric. Nested (write both peptide + protein matrices) proportion **deferred**. |
| **F6** | Low | "warn if proportion high" threshold unspecified. | Lower priority now (decoys dropped at fit → no moderation pollution); keep as a QC proportion metric. |
| **F7** | — | Interim prolfquapp `clean()` pattern-gated-*removes* decoys from quant; deferred design *keeps* them until the model step. | Migration note: flip when implementing; enumerate the `clean()`/`get_subset`/`remove_cont_decoy` sites. |

**All review findings resolved** (F1–F3 in §3; F4/F5 above). Only the nested (peptide + protein
dual-matrix) decoy-proportion case is deferred to a later pass.

## Removed in the prolfquapp target+decoy refactor — reintroduction checklist

The prolfquapp refactor **atomically removes** the REV-column decoy machinery (the root cause is
now fixed by guaranteeing one annotation row per protein ID; decoys are resolved *during*
de-duplication). `pattern_decoys` survives only as a `ProteinAnnotation$new()` arg + the new
`get_rev_pattern()` getter. The inventory below is what was removed and where — use it as the
checklist when reintroducing equivalents as **quant-side** capabilities under the deferred work
above. (Verified by an ecosystem-wide grep on 2026-06-30; line numbers are approximate.)

### API removed from `prolfquapp::ProteinAnnotation` (`R/R6_ProteinAnnotation.R`)

- `annotate_decoys()` (~L302) — the method that created the `REV` column.
- The `REV` column itself (consumed only internally by `clean`/`nr_clean`/`get_summary`).
- The `decoys=` branch of `clean()` (~L379) and `nr_clean()` (~L356).
- Decoy fields in `get_summary()` (~L344-348): `percentOfFalsePositives`, `NrOfProteinsNoDecoys`.

### Kept (NOT removed)

- `pattern_decoys` constructor arg — now drives within-duplicate decoy detection and feeds
  `get_rev_pattern()`.
- All contaminant machinery: `annotate_contaminants()`, `clean(contaminants=)`,
  `nr_clean(contaminants=)`, the `CON` column, `percentOfContaminants`, `remove_cont`.

### Callers updated in prolfquapp

- `R/R6_ProteinDataPrep.R`: `cont_decoy_summary()` (L59-71) loses the two decoy fields, keeps the
  contaminant fields, and gains the construction-time dup/decoy log counts; `remove_cont_decoy()`
  (L78-82) — quant decoy removal becomes **pattern-gated** (per the Q1 decision), contaminant
  removal unchanged.
- `R/cmd_helpers.R`: `cont_decoy_summary()` calls (L345, L570) return fewer fields; IBAQ
  `clean(..., decoys = remove_decoys)` (L458-460) → contaminant-only / pattern-gated.
- `inst/application/CMD_DEA_V2.R`: IBAQ `clean(..., decoys = remove_decoys)` (L222) → same.
- `R/R6_AppConfiguration.R`: the `remove_decoys` boolean (L19) + `remConDec` parsing (L469) —
  `remConDec` still drives contaminant removal; decoy removal is now pattern-gated, so the boolean
  is effectively retired. (Empty-pattern→`NULL` handling already exists at L471.)
- `R/example_deanalyse.R` (L34), `inst/samples/DEAAnalyse.R` (L20): `cont_decoy_summary()` calls —
  still valid, fewer fields.

### Report templates rewritten (drop the "false positives / decoy sequences" sentence)

- `vignettes/Grp2Analysis_V2_R6.Rmd` (L147-148)
- `tests/testthat/Grp2Analysis_V2_R6.Rmd` (L147-148)
- `inst/application/DIANN/_Grp2Analysis_V2.Rmd` (review cites ~L133-137)
- `inst/application/CompoundDiscovery/Grp2Analysis_V2_R6.Rmd` (review cites ~L142-146)
- *Not touched:* generated run-output copies under
  `prolfquasaint/inst/application/SE2_DIANN/.../Inputs_WU_*/Grp2Analysis_V2_R6.Rmd` — historical
  artifacts, not source.

### Tests rewritten

- `tests/testthat/test-R6_ProteinAnnotation.R` (L1-40): the `REV` / `clean(decoys=)` /
  `nr_clean(decoys=)` assertions → replaced by uniqueness + within-duplicate decoy-resolution tests.

### Cross-package impact (small)

- **prophosqua** — no REV-API usage; unaffected.
- **prolfquasaint** — only generated run artifacts reference the decoy summary fields; no source
  change.
- **prolfquappPTMreaders** — passes `pattern_decoys` into `get_annot_from_fasta()` and
  `ProteinAnnotation$new()` (`preprocess_FP_combined_STY.R`, `preprocess_FP_multisite.R`,
  `preprocess_BGS_site.R`); `pattern_decoys` stays valid, but verify behaviour after the
  prolfquapp Stage 3 removes `get_annot_from_fasta()`'s decoy pre-filter (it uses
  `fasta_annot_early`).

### To reintroduce here (the LFQData rev-pattern work)

- **Decoy proportion among quantified entries** (was `percentOfFalsePositives`) — now provided by
  `LFQData$decoy_proportion()` (implemented). Surface it in the prolfquapp summary/QC (edit-list
  item 6). Compute it on the quant data that **still contains decoys** (i.e. any point before the
  fit-drop) — after retiring the `remove_cont_decoy` inner join, that is the peptide LFQData
  throughout, so timing is no longer delicate.
