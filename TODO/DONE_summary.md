# prolfqua — Completed Work Summary

Chronological record of completed development work on the `Modelling2R6` branch.

---

## 2026-03-25 — Item 3: Fix `nrcoef` / `nrcoeff_not_NA` naming inconsistency

From `TODO/TODO_10_refactorings.md` Item 3. Standardized column name spelling from `nrcoef` (one f) to `nrcoeff` (two f's) to match `nrcoeff_not_NA`. Renamed local closure `nrcoeff` → `count_coefficients` to avoid clash with column name.

**Files:** `R/tidyMS_build_model.R`, `R/Model.R`, `R/ModelFirth.R`. All 485 tests pass.

---

## 2026-03-25 — Item 2: Fix cryptic variable names in internal functions

From `TODO/TODO_10_refactorings.md` Item 2. Renamed cryptic 2-letter variables to descriptive names in two files:

**`R/tidyMS_contrasts.R`:**
- `.model_coeff_matrix`: `mm` → `coeff_matrix`, `coefi` → `non_intercept_coeffs`, return field `$mm` → `$coeff_matrix`
- `.get_match_idx`: param `mm` → `coeff_matrix`, `ddd` → `row_name_parts`, `xd` → `factor_match`
- `.coeff_weights_factor_levels`: param `mm` → `coeff_matrix`, `xx` → `weights_by_factor`
- `linfct_from_model`: `cm` → `coeff_result`, `cm_mm` → `sorted_coeff_matrix`

**`R/tidyMS_moderation.R`:**
- `moderated_p_limma`: param `mm` → `contrast_df`, `sv` → `squeezed_var`, `prqt` → `conf_quantile`
- `moderated_p_limma_long`: param `mm` → `contrast_df`, `dfg` → `split_groups`, `xx` → `moderated_results`
- `adjust_p_values`: param `mm` → `contrast_df`, `dfg` → `grouped_df`, `xx` → `adjusted_df`
- `contrasts_fisher_exact`: param `xb` → `fisher_input`, `xx` → `enriched_result`

All 485 tests pass.

---

## 2026-02-19 — Bug Fixes: Hardcoded Imputation, Log Protection, Silent Drops

From code review top-10 items #1, #2, #5.

- **Fixed hardcoded imputation parameters** (`R/LFQDataImp.R`): removed hardcoded overrides, now uses `match.arg()`. Cleaned up dead code, proper return type, fleshed out `LFQDataImp` class.
- **Added log(0)/log(negative) protection** (`R/tidyMS_R6_TransitionCorrelations.R`): `transform_work_intensity()` now warns when log transform encounters zeros (-Inf) or negatives (NaN).
- **Added reporting for silently dropped proteins/contrasts** (`R/tidyMS_R6_Modelling.R`, `R/logistf.R`): `contrasts_linfct()` and `contrasts_linfct_logistf()` now report how many proteins were dropped via `message()`.

---

## 2026-02-20 — Test Coverage & Deprecated Function Migration

Reference: `code_review_report.md`, commit `f1adad05`.

### Added `@examples`

| Class/Function | File | Status |
|----------------|------|--------|
| `LFQDataWriter` — `get_long()`, `get_wide()`, `write_long()`, `write_wide()` | `R/LFQDataWriter.R` | Done |
| `AnalysisConfiguration` — class constructor | `R/AnalysisConfiguration.R` | Done |
| `model_summary()` | `R/tidyMS_R6Model.R` | Done |
| `ionstar_bench_preprocess()` | `R/Benchmark.R` | Done |
| `ms_bench_auc()` | `R/Benchmark.R` | Done |

### Removed `\dontrun{}` wrappers

All `\dontrun{}` removed — code now runs during `R CMD check`:

| File | What was unwrapped |
|------|--------------------|
| `R/LFQDataTransformer.R` | `preprocessCore` quantile normalization (guarded by `if(require(...))`) |
| `R/LFQDataPlotter.R` | `pairs_smooth()` call |
| `R/LFQDataAggregator.R` | `write_plots()` to tempdir |
| `R/tidyMS_missigness.R` | `UpSetR::upset()` interaction plot |
| `R/tidyMS_missigness.R` | `UpSetR::upset()` sample-level plot |

### Cleaned up `tests/testthat/test-plotting_functions.R`

- Removed placeholder test (`"multiplication works"`, `expect_equal(2 * 2, 4)`)
- Replaced `stopifnot("ggplot" %in% class(p))` with `expect_true("ggplot" %in% class(p))` (3 locations)

### Replaced ~92 deprecated dplyr/tidyr/ggplot2 calls across 14 files

| Deprecated | Replacement | Count |
|---|---|---|
| `group_by_at(vars)` | `group_by(across(all_of(vars)))` | 32 |
| `select_at(vars)` | `select(all_of(vars))` | 17 |
| `aes_string(x = var)` | `aes(x = .data[[var]])` | 17 |
| `spread(key, val)` | `pivot_wider(names_from, values_from)` | 6 |
| `one_of(vars)` | `all_of(vars)` | 6 |
| `UQ(sym(x))` | `!!sym(x)` | 5 |
| `unnest_legacy()` | `unnest(cols = c(...))` | 5 |
| `gather(key, val, cols)` | `pivot_longer(cols, names_to, values_to)` | 4 |
| `summarise_at(vars, fn)` | `summarise(across(all_of(vars), fn))` | 4 |
| `mutate_at(var, fn)` | `mutate(across(all_of(var), fn))` | 1 |

Special cases: `aes_string()` with expression strings → `!!rlang::parse_expr()`, nullable colour/text → conditional `aes()`.

### Removed dead code

- `LFQDataWriter` dead `write_plots()` method removed
- Dead `NA <- NA` assignment in `R/LFQDataImp.R` noted

---

## 2026-02-20 — Code Review Items Completed (items 1–7, 9–10)

From the 5-agent code review report (2026-02-19):

| # | Issue | Resolution |
|---|-------|------------|
| 1 | Hardcoded imputation params | `match.arg()` in `impute_with_zcomp()` |
| 2 | log(0)/log(negative) protection | Warnings in `transform_work_intensity()` |
| 3 | 115+ deprecated tidyverse functions | All replaced (commit `f1adad05`) |
| 4 | Missing `@examples` | All added |
| 5 | Silent dropped proteins/contrasts | `message()` in `contrasts_linfct` |
| 6 | MAD division-by-zero guard | Zero-MAD check + warning in `robust_scale()` |
| 7 | Document MCAR assumption | `@note` on `impute_with_zcomp()`, `@section` on `LFQData` |
| 9 | Replace for-loops | Vectorized: `loopOverNested` → `purrr::map()`, `poolvar` → `map_df()`, `percentage_abundance` → grouped mutate, stats crossing → `tidyr::crossing()` |
| 10 | Remove dead code/stubs | Placeholder test, `\dontrun{}`, dead `write_plots()` removed |

Additional: `ContrastsPlotter` code deduplication (`.ma_fig()` helper).

---

## 2026-03 — LFQDataWriter Removal

Commit `3a094331`. `LFQDataWriter` class removed entirely — only used externally by ptm-pipeline, which was updated.

---

## 2026-03 — Limma Backend Added

Commit `cbed2721`. New file `R/ContrastsLimma.R` with:

- `strategy_limma()` — analogous to `strategy_lm()`, returns strategy list with formula, trend, robust
- `build_model_limma()` — fits limma matrix model via `lmFit()` on wide expression matrix
- `ModelLimma` R6 class — wraps limma `MArrayLM` fit, same API as `Model`
- `ContrastsLimma` R6 class — inherits `ContrastsInterface`, uses `contrasts.fit()` + `eBayes()`

User-facing workflow is identical — just swap `strategy_lm` → `strategy_limma` and `Contrasts` → `ContrastsLimma`.

---

## 2026-03 — DEqMS Backend Added

`ContrastsDEqMSFacade` and `ContrastsModeratedDEqMS` implemented. Facade derives count info directly from `LFQData$config$nr_children`. `DEqMS_Moderation.Rmd` vignette added.

---

## 2026-03-12 — Contrast Facades, DEqMS Simplification, Pkgdown Fixes

Commit `dcc7b107`.

### Simplified DEqMS facade API
- Removed `count_df` and `count_column` from `ContrastsDEqMSFacade`
- Facade now derives counts directly from `LFQData` via `lfqdata$config$nr_children`

### Tightened facade input contracts
- Aggregate-only facades (`lm`, `limma`, `deqms`, `lm_missing`) require protein-level `LFQData`
- Nested-data facades (`ropeca`, `lmer`) require peptide-level data
- Added explicit validation errors

### Added `ContrastsLmerFacade`
- Added `method = "lmer"` support in `build_contrast_analysis()`

### Harmonized facade outputs
- Added `facade` column in each facade result
- Kept `modelName` as the underlying engine/result label

### Vignette and test updates
- Added `vignettes/ContrastFacades.Rmd`
- Extended `tests/testthat/test-ContrastsFacades.R`
- Fixed pkgdown metadata/accessibility warnings (site URL, aria-label, alt text)

---

## 2026-03 — Formatting (Air) & Linting Standardization

Air formatter applied. `.lintr` configured with 120-char line length. `object_name_linter` disabled for mixed naming conventions.

---

## 2026-03-18 — Limma Weights Support Added

From `Expose_additional_lmFit_paths.md`. Extended `strategy_limma()` and `build_model_limma()` with `weights` parameter support. Per-sample weights are passed through to `limma::lmFit(weights=)`.

---

## 2026-03-19 — ContrastsProDA Removed, ContrastsRLMFacade Added, Weights in strategy_lm

- `ContrastsProDA` class removed (proDA dependency dropped)
- `ContrastsRLMFacade` added for robust linear model contrasts
- `strategy_lm()` extended with weights support

---

## 2026-03-22 — LM-based Imputation with Covariance Borrowing

From `TODO_new_Imputation.md`. Added a second-pass imputation mechanism inside `build_model()` for proteins whose initial `lm()` fit fails or produces NA coefficients.

### New functions (`R/tidyMS_R6_Modelling.R`)
- `new_lm_imputed()` — S3 constructor wrapping an `lm` object with borrowed covariance; overrides `vcov()`, `sigma()`, `df.residual()` so contrast code works unchanged
- `compute_borrowed_variance()` — computes median sigma or element-wise median vcov from successful fits (two methods: `"sigma"` and `"vcov"`)
- `impute_refit_singular()` — for each failed/singular protein: imputes NAs with LOD, clamps values to `max(value, LOD)`, refits, wraps model with borrowed covariance

### Modified `build_model()` (`R/tidyMS_R6Model.R`)
- New parameters: `impute`, `lod`, `borrow_method` (`"sigma"` / `"vcov"`), `df_method` (`"observed"` / `"borrowed"`)
- When `impute = TRUE`, auto-computes LOD if not supplied, runs `impute_refit_singular()`, appends "Imputed" to modelName

### New facade (`R/ContrastsFacades.R`)
- `ContrastsLMImputeFacade` — same pattern as `ContrastsLMFacade` but passes `impute = TRUE` to `build_model()`

### Tests (`tests/testthat/test-ImputeModel.R`)
- 6 tests, 20 assertions: S3 dispatch, end-to-end model+contrasts, facade comparison, both borrow methods, both df methods

---

## 2026-03-24 — Remove Deprecated `interaction_missing_stats`, Rd Line Width Fixes

Commit `4cb6b74f`.

### Removed `interaction_missing_stats()` (`R/tidyMS_missigness.R`)
- Deprecated function replaced by `summarize_stats()` from `R/tidyMS_stats.R`
- Added `config$isotopeLabel` to `summarize_stats()` grouping for structural parity with the old function
- Migrated all 4 internal callers: `missigness_histogram()`, `missingness_per_condition_cumsum()`, `missingness_per_condition()`, `UpSet_interaction_missing_stats()`
- Deleted function, roxygen block, man page, and NAMESPACE export

### Fixed Rd line width NOTEs
- Wrapped long `hierarchy[["fragment_Id"]]` example lines in `R/AnalysisConfiguration.R` and `R/tidyMS_R6_ConcreteConfigurations.R`
- Eliminated `checking Rd line widths` NOTE from `R CMD check`

---

## 2026-03-24 — Strategy R6 Conversion & Default Weights Enforcement

### Converted all model strategies from plain lists to R6 classes

| Class | File | Wrapper | Fitting function |
|-------|------|---------|-----------------|
| `StrategyLM` | `tidyMS_R6_Modelling.R` | `strategy_lm()` | `lm()` |
| `StrategyRLM` | `tidyMS_R6_Modelling.R` | `strategy_rlm()` | `MASS::rlm()` |
| `StrategyLmer` | `tidyMS_R6_Modelling.R` | `strategy_lmer()` | `lmerTest::lmer()` |
| `StrategyLogistf` | `logistf.R` | `strategy_logistf()` | `logistf::logistf()` |

Each class exposes: `formula`, `model_name`, `report_columns`, `is_mixed`, `anova_df` (fields) and `model_fun()`, `isSingular()`, `contrast_fun()`, `df_residual()`, `sigma()` (methods). Wrapper functions preserved for backward compatibility — all existing call sites unchanged.

### Dropped `strategy_glm` (dead code)
- Zero production callers, zero test callers — only referenced in its own `@examples`

### Enforced `nr_children` as default fitting weights

**Key distinction clarified:** `nr_children` is **sample-wise** (per protein×sample count of child features after aggregation), not experiment-wide.

- 5 aggregated facades now accept `weights = lfqdata$config$nr_children` by default:
  `ContrastsLimmaFacade`, `ContrastsLMFacade`, `ContrastsLMMissingFacade`, `ContrastsLMImputeFacade`, `ContrastsDEqMSFacade`
- Users pass `weights = NULL` to disable
- `StrategyLM` holds `weights` field; passed to `lm(..., weights=)` via `bquote`
- `build_model_limma()` weight handling fixed: columns in `value_vars()` (like `nr_children`) pivoted to protein×sample matrix via `to_wide()`; custom per-sample columns use vector extraction

**DEqMS uses `nr_children` differently:** aggregates via `max()` per protein (experiment-wide) for count-dependent variance moderation — separate from fitting weights.

### Files modified

| File | Changes |
|------|---------|
| `R/tidyMS_R6_Modelling.R` | `StrategyLM`, `StrategyRLM`, `StrategyLmer` R6 classes; dropped `strategy_glm` |
| `R/logistf.R` | `StrategyLogistf` R6 class |
| `R/ContrastsFacades.R` | `weights` parameter added to 5 aggregated facades |
| `R/ContrastsLimma.R` | `build_model_limma` weight pivot for protein×sample matrix |
| `CLAUDE.md` | Updated strategy pattern and weights documentation |

---

## 2026-03-24 — File Reorganization: Renames, Splits, Dead Code Removal

### Removed dead files and moved functions to their consumers

Commits `edc5dee3`, `c305a1ce`, `a1d79433`.

| Deleted file | Reason | Where functions went |
|---|---|---|
| `R/tidyMS_R6_ConcreteConfigurations.R` | `create_config_MQ_peptide()` dead (zero callers) | Deleted entirely |
| `R/tidyMS_MQ_workflow.R` | `filter_proteins_by_peptide_count()`, `filter_difference()` only called from `LFQData` | Moved to `R/LFQData.R` |
| `R/tidyMS_R6_TransitionCorrelations.R` | Mixed bag of unrelated functions | Split by theme (see below) |

### Split `tidyMS_R6_TransitionCorrelations.R` by theme (commit `c305a1ce`)

| Destination | Functions moved |
|---|---|
| `R/LFQData.R` | `remove_small_intensities`, `nr_B_in_A`, `.nr_B_in_A`, `.make_name_AinB` |
| `R/LFQDataTransformer.R` | `transform_work_intensity`, `response_matrix_as_tibble`, `.get_robscales`, `get_robscales`, `robust_scale`, `apply_to_response_matrix`, `scale_with_subset`, `center_to_reference`, `center_to_reference_cfg` |
| `R/tidyMS_aggregation.R` | `.ExtractMatrix`, `nr_obs_sample`, `.rankProteinPrecursors`, `rank_peptide_by_intensity` |
| `R/tidyMS_reshaping.R` (new) | `tidy_to_wide`, `tidy_to_wide_config` |

Also removed dead `nr_obs_experiment()` (zero callers).

### Rename `tidyMS_missigness_V2.R` → `tidyMS_missigness_imputation.R` (commit `a1d79433`)

Name now reflects actual content (imputation helpers, not "V2" of missingness).

### Rename `tidyMS_R6Model.R` → `tidyMS_build_model.R` and split `tidyMS_R6_Modelling.R` (commit `d04e4775`)

`tidyMS_R6Model.R` had no R6 classes — renamed to `tidyMS_build_model.R` to match its content (`build_model`, `build_model_impute`, `LR_test`, `model_summary`).

`tidyMS_R6_Modelling.R` (1634 lines, 5 unrelated concerns) split into 4 files:

| File | Lines | Contents |
|---|---|---|
| `R/tidyMS_R6_Modelling.R` | 396 | Strategy R6 classes (`StrategyLM`, `StrategyRLM`, `StrategyLmer`, `AnovaExtractor`) + wrappers |
| `R/tidyMS_build_model.R` | 605 | Model fitting internals merged in: `model_analyse`, `isSingular_lm`, `get_complete_model_fit`, `new_lm_imputed` + S3 methods, `compute_borrowed_variance`, `impute_refit_singular`, `plot_lmer_peptide_predictions` |
| `R/tidyMS_contrasts.R` (new) | 548 | Contrast linear function machinery: `linfct_*` family, `my_contrast`, `my_contrast_V2`, `my_contest`, `contrasts_linfct`, `pivot_model_contrasts_2_Wide` |
| `R/tidyMS_moderation.R` (new) | 318 | Moderation + ROPECA + Fisher: `moderated_p_limma*`, `adjust_p_values`, `get_p_values_pbeta`, `summary_ROPECA_median_p.scaled`, `contrasts_fisher_exact` |

No API changes — all exports preserved, all 485 tests pass.
