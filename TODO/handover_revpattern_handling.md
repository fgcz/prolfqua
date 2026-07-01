# Handover — decoy / contaminant handling (rev-pattern work)

**Date:** 2026-07-01
**Author:** Claude (agent session)
**For:** review of the changes across `prolfqua` and `prolfquapp`.
**Companion design doc:** `prolfqua/TODO/TODO_revpattern_handling.md` (the *why* + locked
decisions + adversarial review). This file is the *what changed / how to review / how to verify*.

---

## 1. TL;DR

The original bug (workunit **WU347806**) was a target+decoy FASTA producing duplicate `protein_Id`s,
which crashed `column_to_rownames()` during SummarizedExperiment export. That root cause was fixed
earlier by making `ProteinAnnotation` guarantee one row per protein id (`prolfquapp c9d6ddb`).

This follow-up work then:

1. **Promoted decoy/contaminant detection into `prolfqua` core** (one shared, database-agnostic
   detector) and gave `LFQData` a **targets-only model fit** so decoys can never pollute the shared
   variance pool.
2. **Rebuilt prolfquapp onto a single explicit quant filtering path** — removing a hidden/redundant
   filter (an annotation-alignment inner join that silently dropped decoys) and settling on:
   - **decoys** — kept through normalization, dropped **only at the model fit**, preserved in the
     raw/abundance data for export (they carry NA statistics → invisible in results/volcano/ORA);
   - **contaminants** — **kept everywhere and labelled**, never removed.
3. **Added interactive report plots (ggiraph)** — a contaminant-marked volcano and per-sample
   intensity-density curves.

**Net behaviour change to flag:** contaminants now **appear in the results table and volcano
(flagged), instead of being removed** — a deliberate flip of the previous default. Decoys now appear
in the abundance matrix / SummarizedExperiment with NA statistics (previously removed upstream).

---

## 2. Design decisions (locked with you)

| Topic | Decision |
|---|---|
| Decoy detection | shared `prolfqua::is_decoy()` (config pattern ∪ built-in defaults; empty/`NULL`/`"a^"` → defaults only, never `grepl("", x)`) |
| Where decoys are dropped | **only at the model fit** (targets-only), so the limma prior / DEqMS variance-count trend never sees them (**F1**) |
| Decoys on export | kept in the abundance matrix / SE with **NA** diff/p/FDR → absent from volcano / ORA / GSEA / significant set automatically (**F2**) |
| Contaminants | **kept + labelled, never removed** (real proteins); flagged via the annotation `CON` flag, which rides the export `right_join` onto results |
| `LFQData$remove_contaminants()` | **deleted** (removal is not part of the pipeline) |
| Filtering paths | **one** — the `get_subset(clean())` inner join in `ProteinDataPrep`/IBAQ is retired; the export `right_join` is the sole annotation↔quant join |
| Normalization | runs on the full data (decoys + contaminants included); the ~1% shift is accepted (**F3**) |

---

## 3. Commits to review

### prolfqua (branch `main`)

| Commit | Summary |
|---|---|
| `f0f89aec` | Design note (initial) |
| `706e08af` | Shared detectors: `is_decoy` / `is_contaminant` / `effective_*_pattern` (`R/utilities.R`) + tests |
| `4166581f` | `AnalysisConfiguration` `pattern_decoys`/`pattern_contaminants` slots; `LFQData$remove_decoys`/`decoy_proportion`/`contaminant_proportion`; targets-only fit in `build_contrast_analysis` + tests |
| `b046e9c3`, `1df4158f`, `8e3361d1`, `0b13710d` | design-doc updates (status, single-path investigation, keep+label) |
| `b85b535e` | **Delete `LFQData$remove_contaminants()`** (keep+label); test updated |

### prolfquapp (branch `master`)

| Commit | Summary |
|---|---|
| `c9d6ddb` | (root-cause fix) `ProteinAnnotation` = single unique-id + decoy authority |
| `e3a2bcd` | Consolidate: `.detect_decoy_ids` delegates to `prolfqua::is_decoy` (removes duplicate) |
| `fdbf171` | **Single filtering path**: Gap B (stamp patterns on modelling config), Gap A (drop decoys in `DEAnalyse$build_facade`), retire `remove_cont_decoy` filter, decoy QC in `cont_decoy_summary`; new `test-single-filtering-path.R` |
| `f46eeeb` | Align IBAQ to the single path (no `get_subset(clean())` pre-filter) |
| `bd2e495` | Interactive ggiraph report plots (contaminant-marked volcano + density) + `test-interactive-report-plots.R` |

---

## 4. Files to review (by concern)

**prolfqua core**
- `R/utilities.R` — `is_decoy` / `is_contaminant` / `effective_decoy_pattern` / `effective_contaminant_pattern` (+ `.default_*_prefixes`, `.effective_prefix_pattern`).
- `R/AnalysisConfiguration.R` — `pattern_decoys` / `pattern_contaminants` fields (default `NULL`; round-trip through `R6_extract_values`).
- `R/LFQData.R` — `remove_decoys()`, `decoy_proportion()`, `contaminant_proportion()`, private `.prefix_proportion()`.
- `R/build_contrast_analysis.R` — gated targets-only drop (opt-in via `pattern_decoys`).

**prolfquapp filtering path**
- `R/R6_ProteinAnnotation.R` — `.detect_decoy_ids` now delegates to `prolfqua::is_decoy`.
- `R/R6_ProteinDataPrep.R` — `initialize` stamps patterns; `cont_decoy_summary()` adds `percentOfDecoys`; **`remove_cont_decoy()` removed**.
- `R/R6_DEAnalyse.R` — `initialize` stamps patterns onto the modelling config; **`build_facade` drops decoys before the fit** (covers nested subclass via inheritance + SAINT branch).
- `R/cmd_helpers.R` — two `remove_cont_decoy()` calls removed; IBAQ uses the full LFQData (a copy).

**prolfquapp reporting (ggiraph)**
- `R/interactive_report_plots.R` — `volcano_ggiraph()`, `intensity_density_ggiraph()`.
- `vignettes/Grp2Analysis_V2_R6.Rmd` — density chunk + volcano chunk + caption.
- `DESCRIPTION` — `ggiraph` added to Suggests.

**Tests added**
- prolfqua: `test-decoy-contaminant-detect.R`, `test-decoy-contaminant-lfqdata.R`, `test-build-contrast-decoy-drop.R`, config round-trip in `test-tidyconfig_functions.R`.
- prolfquapp: `test-single-filtering-path.R`, `test-interactive-report-plots.R`.

---

## 5. How to verify

```bash
# prolfqua unit tests (expect 1047 pass, 0 fail)
cd prolfqua && Rscript -e "devtools::load_all('.'); testthat::test_dir('tests/testthat')"

# install prolfqua so prolfquapp sees the updated core
cd prolfqua && make install

# prolfquapp unit tests (expect 257 pass, 0 fail)
cd prolfquapp && Rscript -e "devtools::load_all('.'); testthat::test_dir('tests/testthat')"
```

**WU347806 end-to-end** (real DIA-NN data): re-verified — SE builds, 4029 unique proteins, RDS
written. The run script is in the session scratchpad (`run_se_validation.R`); it does
`get_config → sync_opt_config → run_dea → make_SummarizedExperiment`.

**Report render** (the two ggiraph chunks): verified rendering to self-contained HTML with embedded
girafe widgets via a minimal driver Rmd built from a synthetic decoy+contaminant fixture
(`scratchpad/render_ggiraph_chunks.R`). To eyeball the plots, open that HTML, or run a full DEA and
open the generated `Grp2Analysis_V2_R6.html`.

---

## 6. What the fixture proves (`test-single-filtering-path.R`)

Using `add_RevCon()` to inject `REV_` (~10%) and `zz` (~5%) onto simulated peptide data:
- `cont_decoy_summary()` reports `percentOfDecoys > 0` and **does not filter**;
- decoys **and** contaminants survive aggregation + normalization;
- after the fit: **no decoy in the contrasts**, contaminants **present** in the contrasts, and
  decoys still present in `lfq_data_raw` (for export);
- `pattern_decoys = NULL` (opt-out) → decoy machinery off, decoys **not** dropped
  at the fit; `pattern_decoys = ""` opts **in** (defaults only) and decoys **are**
  dropped — same semantics as `prolfqua::build_contrast_analysis` (the app maps an
  empty REVpattern to NULL upstream).

---

## 7. Open / deferred (not done — your call)

- **Nested (peptide + protein) decoy-proportion** — deferred (design doc §"Proportion for nested").
- **`remove_cont` config flag** — now vestigial for the default path; kept in case a future opt-in
  "biologist removes contaminants" mode is wanted (additive).
- **FDR-validation volcano** (decoys should be non-significant) — not available, because decoys are
  not modelled. Would need a *separate* fit that does not share the target variance pool.
- **Blocker C (DIA-NN strips `zz|` from the quant id)** — sidestepped, not fixed: contaminant
  *labelling* uses the annotation `CON` flag (works for DIA-NN), so we never needed to detect
  contaminants on the DIA-NN quant id. If a future feature needs quant-side contaminant detection
  for DIA-NN, revisit `preprocess_DIANN.R:91`.

## 8. Watch-items / risks

1. **Pattern must reach the modelling config.** The targets-only drop is gated on
   `config$pattern_decoys`; it is stamped in `DEAnalyse$initialize` (robust — after aggregation).
   If a new code path constructs a facade without going through `DEAnalyse`, decoys could re-enter
   the fit. All current paths (`build_deanalyse`, nested, SAINT) are covered.
2. **Contaminants now in the fit + normalization.** High-abundance contaminants (keratin, etc.)
   participate in normalization and get real contrasts. This is the intended keep+label behaviour,
   but it changes results vs. the old "remove contaminants" default — worth a sanity check on a
   contaminant-heavy dataset.
3. **`ggiraph` is a Suggests dep.** The report helpers `stop()` with a clear message if it is
   missing. Ensure report-rendering environments have `ggiraph` installed.
4. **SE now has more rows** (decoys with NA stats; contaminants with real stats). Any consumer
   assuming `nrow(SE) == nrow(contrasts)` needs auditing; the shipped consumers were checked (ORA
   background uses the decoy-free annotation; the Quarto SE report `na.omit()`s).
