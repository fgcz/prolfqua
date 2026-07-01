# Rev/decoy pattern handling in prolfqua core (deferred)

**Date:** 2026-06-30
**Status:** Deferred / future. **Do not implement yet** — `prolfqua` core stays untouched for now.
**Driver:** the prolfquapp redesign in
`prolfquapp/TODO/TODO_protect_prolfquapp_against_target_decoy.md` (handling target+decoy FASTAs,
guaranteeing unique protein IDs). This note records the `prolfqua`-side follow-up so it is not
lost.

## Background

`prolfqua`'s `LFQData` / `AnalysisConfiguration` currently has **no** decoy/rev concept at all
(the only `rev` is "reverse hierarchy order", unrelated). All decoy knowledge lives in
prolfquapp's `ProteinAnnotation`.

In the prolfquapp redesign, `ProteinAnnotation`:
- is always **decoy-free** (decoys never belong in a protein annotation),
- **detects** the decoy/rev pattern (configured pattern, else a built-in default prefix set),
- **exposes** the detected pattern via a getter (e.g. `get_rev_pattern()`).

## Design (LOCKED 2026-07-01) — quant-side rev/contaminant handling

**Status of this section:** design agreed; **pending one adversarial review** before implementation.
Still deferred (prolfqua core untouched until we pick this up).

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

### 4. Contaminants — annotate + pattern-gated, config-controlled removal (asymmetric to decoys)

- `LFQData` flags `is_contaminant` **only when `pattern_contaminants` is configured** — no pattern
  → contaminants can't be identified → none flagged, none removed. (Matches today's behaviour:
  the `a^` no-op pattern flags no `CON`, so nothing is removed.)
- Contaminants are **real** proteins (keratin, trypsin, BSA), can be high-abundance and distort
  results → a genuine scientific keep/remove choice (the bioinformatics-vs-biologist split).
- When a pattern identifies them, `remove_cont` (**default = remove**) controls the analysis:
  default removes (biologist); flip to keep+flag (bioinformatics).
- **QC:** report contaminant proportion.

### 5. What moves out of `ProteinAnnotation` when we implement this

The interim prolfquapp state routes quant contaminant removal through
`ProteinAnnotation$clean(contaminants=)` + `get_subset`, and pattern-gated decoy removal through
`clean()`. This deferred work **relocates both to `LFQData`** (annotate + optional filter), leaving
`ProteinAnnotation` using the pattern **only** for dedup.

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
| Ever drop contaminants? | only if `pattern_contaminants` set; then **default remove** (`remove_cont`), flip to keep+flag |
| Decoy-proportion QC | `nr_decoys / nr distinct hierarchy keys` of the LFQData (per-peptide or per-protein by analysis level) |
| Proportion for nested (peptide+protein) analysis | **deferred** |

## Why deferred

- Touching `prolfqua` core ripples to all dependents (prolfquapp, prophosqua, …); not warranted to
  fix the immediate prolfquapp bug.
- prolfquapp owns detection now and exposes the pattern; promoting it to core is the clean
  follow-up now that the prolfquapp design has settled.

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

- **Decoy proportion among quantified entries** (was `percentOfFalsePositives`) — recompute on the
  quant side via `LFQData` + `ProteinAnnotation$get_rev_pattern()`, *before* the pattern-gated
  removal, and surface it in QC.
