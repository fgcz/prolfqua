# prolfqua — Completed Work Summary

Chronological record of completed development work on the `Modelling2R6` branch.

---

## 2026-06-23 — TODO Cleanup

Archived `TODO_LFQData_access_patterns.md` after confirming the note itself
marks all LFQData API, standalone-function refactoring, Transformer refactoring,
and private `data`/`config` encapsulation phases as DONE/RESOLVED. The root
`TODO.md` now keeps only the remaining broader encapsulation follow-ups.

---

## 2026-05-28 — TODO Cleanup

Archived completed TODO notes after checking their planned fixes against the current codebase:

- version bump to 1.6.2
- GitHub tarball symlink fix
- limma sample/annotation alignment fix
- Windows limpa vignette compatibility fix
- plotting label and volcano plot fixes
- contrast table score, rank/ORA, and model-color TODOs
- limpa same-facade / lmer nested bug follow-up

Kept active TODOs in the top-level `TODO/` directory when the underlying work still appears open or partially open.

---

## 2026-04-10 — LFQData API & Config Decoupling (Phase 4 complete)

Major refactoring to decouple standalone functions from AnalysisConfiguration. Functions no longer receive `config` — they get either individual column name arguments (≤5 fields) or an `lfqdata` object (>5 fields).

### Phase 1+3: LFQData API methods + decorator updates

Added 9 accessor methods to LFQData: `hierarchy_keys()`, `relevant_hierarchy_keys()`, `factor_keys()`, `relevant_factor_keys()`, `sample_name()`, `file_name()`, `nr_children_col()`, `isotope_label()`, `get_data()`.

Updated LFQDataPlotter, LFQDataStats, LFQDataAggregator to use API methods instead of `self$lfq$config$...` for direct config field reads.

### Phase 4: Standalone function signature refactoring (8 batches)

**Batch 0+1 — Trivial 1-field functions:**
- `remove_small_intensities(pdata, config)` → `(pdata, response, threshold)`
- `sample_subset(size, pdata, config)` → `(size, pdata, hierarchy_keys_depth)`

**Batch 2 — Hierarchy/factor summaries:**
- `table_factors(pdata, config)` → `(pdata, file_name, sample_name, factor_keys)`
- `table_factors_size(pdata, config)` → `(pdata, file_name, sample_name, factor_keys, factor_keys_depth)`
- `hierarchy_counts(pdata, config)` → `(pdata, hierarchy_keys, isotope_label)`
- `summarize_hierarchy(pdata, config)` → `(pdata, hierarchy_keys, isotope_label, hierarchy, factors)`

**Batch 3+4 — Stat plots + intensity distribution:**
- `plot_stat_density/density_median/violin/violin_median` — `config` → `factor_key` or `factor_keys_depth`
- `plot_stdv_vs_mean` — `config` → `factor_keys_depth`
- `plot_intensity_distribution_violin/density` — `config` → `(sample_name, response, is_transformed)`

**Batch 5+7 — Stats core + missingness (merged, accept lfqdata):**
- `complete_cases(pdata, config)` → `(lfqdata)` + internal `.complete_cases_impl` for bootstrap
- `summarize_stats(pdata, config, ...)` → `(lfqdata, factor_key)` — dropped `.completed` flag
- `summarize_stats_all/factors` → `(lfqdata)`
- `summarize_stats_quantiles(stats_res, config)` → `(stats_res, factor_keys_depth)`
- `encode_bin_resp(pdata, config)` → `(lfqdata)` — `config$bin_resp` side-effect moved to callers
- `missigness_histogram`, `upset_missing_stats`, `missingness_per_condition`, `missingness_per_condition_cumsum`, `upset_interaction_missing_stats` — all `(data, config)` → `(lfqdata)`
- `lfq_power_t_test_quantiles(pdata, config)` → `(lfqdata)`

**Batch 6 — Heatmap/PCA plots (pre-computed wide data):**
- `plot_heatmap_cor/heatmap/raster/na_heatmap` — `(data, config)` → `(matrix, annotation, factor_keys, sample_name)`
- `plot_pca` — `(data, config)` → `(matrix, annotation, sample_name, factor_keys)`
- `plot_sample_correlation` — `(pdata, config)` → `(matrix)`
- LFQDataPlotter now calls `self$lfq$to_wide(as.matrix = TRUE)` and passes decomposed parts

**Batch 8 — Aggregation + nr_obs + boxplots:**
- `estimate_intensity(data, config, .func)` → `(lfqdata, .func)` — `.func` callback interface changed
- `medpolish_estimate_dfconfig/rlm_estimate_dfconfig` → individual args `(pdata, response, hierarchy_keys, hierarchy_keys_depth, sample_name)`
- `aggregate_intensity_topN` → renamed to `aggregate_intensity_top_n`, accepts `(ranked_data, lfqdata, .func, N)`
- `nr_obs_sample/nr_obs_experiment` → individual args
- `rank_peptide_by_intensity` → `(pdata, response, hierarchy_keys)`
- `plot_estimate` → `(lfqdata, lfqdata_agg)`
- `plot_hierarchies_line/line_df` → `(pdata, lfqdata)`
- `plot_hierarchies_boxplot/boxplot_df` → `(pdata, lfqdata)`
- `.reestablish_condition`, `.add_nr_children`, `plot_hierarchies_add_quantline` → individual args

**Cross-package (prolfquapp):**
- `R6_ProteinAnnotation.R`: `nr_obs_experiment` call updated to individual args
- `R6_DEAReportGenerator.R`: `upset_interaction_missing_stats` call updated to `(lfqdata)`
- Vignette `Grp2Analysis_V2_R6.Rmd` + `doc/` + `inst/doc/` copies updated

### Design decisions recorded

- Functions needing >5 config fields get `lfqdata`; otherwise data frame + column names
- `complete_cases` runs at LFQData lifecycle only — dropped `.completed` flag from all functions
- Heatmap/PCA functions receive pre-computed wide data from decorator
- `get_data()` → `data_long()` rename deferred to avoid churn

### Verification

- prolfqua: `make check` (full, with vignettes) — 0 errors, 0 warnings
- prolfquapp: `devtools::check(vignettes = FALSE)` — 0 errors, 0 warnings

### Remaining

- Phase 5: Make `data`/`config` private with active bindings
- Deferred: `get_data()` → `data_long()` + add `data_wide()`

---

## 2026-04-12 — Phase 2: Refactor Transformer to Return New LFQData

Each transform method now creates a new LFQData instance instead of mutating a clone's data and config in place. Constructor stores reference (not clone). User-facing API unchanged: `$get_Transformer()$log2()$robscale()$lfq`.

### Methods refactored

| Method | Change |
|--------|--------|
| `initialize` | Reference instead of deep clone |
| `log2` | Clone config → `transform_work_intensity` → new LFQData |
| `intensity_array` | Same pattern as log2 |
| `intensity_matrix` | Clone config → `apply_to_response_matrix` → new LFQData |
| `robscale` / `robscale_subset` | Clone config → `scale_with_subset` → rename → new LFQData |
| `center_to_reference` | Switch from `copy = FALSE` to `copy = TRUE` → store result |
| `get_scales` | Read-only, no changes |

### Verification

- `make check-fast` — 0 errors, 0 warnings, 0 notes
- No caller changes needed — API identical

---

## 2026-04-07/08 — Complete camelCase to snake_case migration (Phases 1–5)

Completed the full snake_case migration across prolfqua and downstream packages (prolfquapp, prolfquasaint, prophosqua). 9 commits over 2 days.

### Phase 1: Local variables (done prior)

All camelCase local variables renamed across ~20 R files. Inner helper functions (`getCoeffs`, `getSampleSize`, `getAST`) also renamed.

### Phase 2: Function parameters

Renamed `modelName` → `model_name` parameter across 19 functions/constructors (16 R files). Renamed `modelDF` → `model_df` parameter in `compute_borrowed_variance()` and `impute_refit_singular()`. Also renamed `preserveMean`, `sampleName`, `proteinName`, `factorDepth`, `nrNA`, `intesityNewName`, `modelName_Int`, `modelWithInteractionsContrasts` (done in prior commit).

### Phase 3: Exported function names

| Old | New |
|-----|-----|
| `scriptCopyHelperVec` | `script_copy_helper_vec` |
| `pivot_model_contrasts_2_Wide` | `pivot_model_contrasts_to_wide` |
| `isSingular_lm` | `is_singular_lm` |
| `sim_lfq_data_2Factor_config` | `sim_lfq_data_2factor_config` |
| `UpSet_interaction_missing_stats` | `upset_interaction_missing_stats` |
| `get_UniprotID_from_fasta_header` | `get_uniprot_id_from_fasta_header` |

All internal callers, tests, vignettes, downstream packages updated. NAMESPACE regenerated.

### Phase 4a: Delete deprecated R6 aliases

Removed `hierarchyKeys()` and `hkeysDepth()` from AnalysisConfiguration. Fixed 2 internal callers in LFQData.R and 2 downstream callers in prolfquapp.

### Phase 4b: R6 method renames

| Old | New | Class |
|-----|-----|-------|
| `subject_Id()` | `subject_id()` | LFQData |
| `omit_NA()` | `omit_na()` | LFQData |
| `NA_heatmap()` | `na_heatmap()` | LFQDataPlotter |
| `get_LOD()` | `get_lod()` | MissingHelpers |
| `plot_NA_heatmap()` | `plot_na_heatmap()` | standalone function |

### Additional parameter/local renames

- `FCthreshold`/`FDRthreshold` → `fc_threshold`/`fdr_threshold` (8 `get_Plotter` methods)
- `isotopeLabel` param → `isotope_label`, `sampleName` param → `sample_name`
- `remove_NA_rows` → `remove_na_rows`, `UpSet_missing_stats` → `upset_missing_stats`
- ~20 additional local variables: `relativeRisk`, `odsRatio`, `apply_fischer`, `hierarchy_ID`, `summaryColumn`, `rankColumn`, `PCx`/`PCy`, `forPairs`, `levelA`/`levelB`, `notNA`, `resData`, `resMat`, `fileName` param, `maxNrOfSignificantText`, etc.

### Phase 5: R6 field names

| Old field | New field | Classes affected |
|-----------|-----------|-----------------|
| `contrastDF` | `contrast_df` | ContrastsPlotter |
| `modelDF` | `model_df` | Model (+ all callers) |
| `modelName` | `model_name` | 12 classes (Contrasts, ContrastsLimma, Model, ModelFirth, ContrastsModerated, ContrastsDEqMS, ContrastsMissing, ContrastsROPECA, ContrastsFirth, ContrastsTable, ContrastsPlotter, ModelLimma) |
| `subject_Id` | `subject_id` | Same 12 classes + function params in tidyMS_build_model.R, tidyMS_contrasts.R, tidyMS_moderation.R, logistf.R |

**Critical:** `"modelName"` string literals (column names in output data frames) preserved unchanged.

### Bug fix

Fixed `rlm_estimate` roxygen example: built-in datasets use column `sampleName` (data frame column, not renamed), but examples incorrectly passed `"sample_name"`.

### Verification

- prolfqua: `make check-fast` → 0 errors, 0 warnings
- prolfquapp: `make check-fast` → 0 errors, 0 warnings
- prolfquasaint, prophosqua: all renamed callers updated, no stale references

### Remaining camelCase (deliberately kept)

- **R6 class names** — PascalCase by convention (`AnalysisConfiguration`, `ContrastsLimma`, etc.)
- **Factory methods** — reference class names (`get_Plotter`, `get_Transformer`, `get_Aggregator`, etc.)
- **R stats conventions** — `p.value`, `std.error`, `conf.low`, `p.adjust`, `sig.level`, `sample.var`, etc.
- **squeezeVarRob.R** — ported from limma, matches upstream naming
- **Simulation params/constants** — `Nprot`, `N`, `FC`, `PEPTIDE`, `TWO`
- **Math/stats variables** — `Sigma.hat`, `X`, `M`, `SS`
- **Data frame column name strings** — `"modelName"`, `"sampleName"`, `"isotopeLabel"`, etc.

### Still outstanding (not part of snake_case migration)

- Phase 6 (encapsulation): Add accessors for `hierarchy`, `factors` on AnalysisConfiguration; encapsulate `$data`/`$config` on LFQData — deferred
- Phase 7 (`(data, config)` → `(lfqdata)` signatures): Large refactor — deferred
- Phase 8 cross-package: 2 `AnalysisTableAnnotation$new()` references in prolfquasaint (DIANN_SE.R:54, CreateSaintExpress_Report.R:63) — should be `AnalysisConfiguration$new()`

---

## 2026-03-31 — Rename AnalysisConfiguration fields from camelCase to snake_case

Normalized all 9 remaining camelCase/mixed-case fields in `AnalysisConfiguration` to snake_case: `fileName` → `file_name`, `sampleName` → `sample_name`, `normValue` → `norm_value`, `isotopeLabel` → `isotope_label`, `workIntensity` → `work_intensity`, `factorDepth` → `factor_depth`, `hierarchyDepth` → `hierarchy_depth`, `ident_qValue` → `ident_q_value`, `ident_Score` → `ident_score`.

### Active binding aliases for backward compatibility

Added 9 active bindings so old camelCase names still work (read + write) for downstream packages (prolfquapp, prophosqua, etc.). No deprecation warnings yet.

### R6_extract_values() — serialization fix

Added filtering via `.__enclos_env__$.__active__` to exclude active binding aliases from serialized output, preventing duplicate keys.

### Constructor simplified

Dropped `analysisTableAnnotation` / `analysisParameter` parameters. Constructor is now `initialize = function() {}`.

### Dropped `table` / `parameter` active bindings

Removed the deprecated `config$table` and `config$parameter` aliases that returned `self`.

### Files modified

| Area | Files | Changes |
|------|-------|---------|
| Core | `R/AnalysisConfiguration.R` | Field renames, 9 active bindings, `R6_extract_values()` fix, `make_reduced_hierarchy_config()` update |
| R/ files | ~20 files across R/ | ~242 references migrated from camelCase to snake_case |
| Tests | `tests/testthat/test-tidyconfig_functions.R`, `test-Contrasts.R`, `test-plotting_functions.R` | Config field references updated |
| Vignettes | `CreatingConfigurations.Rmd`, `SimulateData.Rmd`, `ContrastFacade2Factor.Rmd`, `LimmaBackend.Rmd` | Config field references updated |
| Other | `README.md`, `CLAUDE.md`, `data-raw/fix_deprecated_config.R`, `inst/issue71/prolfqua_contrast_testCode.R` | Config field references updated |

### Downstream impact

- Active binding aliases keep downstream packages working — no immediate breakage for field renames
- `$new(atable)` constructor pattern broken (~17 calls in prolfquapp/prolfquappPTMreaders/prolfquasaint)
- `config$table$X` / `config$parameter$X` now errors — fix: `config$X`

All 618 tests pass (0 failures, 42 pre-existing warnings).

---

## 2026-03-30 — Limpa Integration: AggregateLimpa, build_model_limpa, ContrastsLimpaFacade

Full integration of the limpa package (Li & Smyth) into the prolfqua modelling pipeline. Limpa provides probabilistic protein quantification from precursor-level data with missing value handling via a Detection Probability Curve (DPC).

### AggregateLimpa (`R/LFQDataAggregator.R`)

New R6 aggregator class following the existing `AggregateMedpolish`/`AggregateRlm`/`AggregateTopN` pattern:

- **Aggregation mode** (`impute_only = FALSE`): wraps `limpa::dpc()` + `limpa::dpcQuant()` to aggregate peptides → proteins. Output is a protein-level `LFQData` with no NAs, plus SE and observation count columns.
- **Impute-only mode** (`impute_only = TRUE`): wraps `limpa::dpc()` + `limpa::dpcQuantByRow()` to fill missing values at the same hierarchy level (e.g. peptide-level imputation for PTM analysis). Hierarchy is preserved.

Helper `.elist_to_lfqdata()` converts limpa's EList output back to long-format LFQData, parsing `~lfq~`-separated rownames back to hierarchy columns and attaching SE + n_obs as extra columns via `config$opt_se` and `config$nr_children`.

### AnalysisConfiguration: `opt_se` field (`R/AnalysisConfiguration.R`)

Added `opt_se` field for standard error column names, following the `opt_rt`/`opt_mz` pattern. Three one-line edits: field declaration, `initialize()` fields list, `value_vars()` return.

### StrategyLimpa + build_model_limpa (`R/ContrastsLimpa.R`)

- `StrategyLimpa` R6 class with fields for formula, model_name, trend, robust, plot, span
- `strategy_limpa()` factory function
- `build_model_limpa()`: extracts SE matrix → `log(SE + 1e-6)` as vooma predictor, n_obs matrix → imputed flag `(n_obs == 0)`, calls `limpa::voomaLmFitWithImputation()` for bivariate vooma precision weighting, returns `ModelLimma`

Since the output is standard `MArrayLM`, `ContrastsLimma` is reused as-is.

### ContrastsLimpaFacade (`R/ContrastsFacades.R`)

New facade class wiring `strategy_limpa` → `build_model_limpa` → `ContrastsLimma`. Validates that `config$opt_se` is set (ensuring `AggregateLimpa` was used). Registered as `method = "limpa"` in `build_contrast_analysis()`.

### Extract `.lfqdata_to_elist()` and `.resolve_weights()` helpers (`R/ContrastsLimma.R`)

Deduplicated the common preamble across all 5 `build_model_limma*` functions:

- `.lfqdata_to_elist(lfqdata, formula)` — `to_wide()` → expr_matrix + annotation + rowdata, subject_Id/isotopeLabel dedup, design matrix, EList construction, dummy lm for linfct extraction
- `.resolve_weights(lfqdata, strategy, annotation)` — resolves `strategy$weights` from column name (per-sample or per-protein×sample) or matrix to the right format for `limma::lmFit()`

Refactored: `build_model_limma`, `build_model_limma_impute`, `build_model_limma_voom`, `build_model_limma_voom_impute`, `build_model_limpa`.

### Bug fix: `.elist_to_lfqdata` duplicate isotopeLabel columns

The annotation join in `.elist_to_lfqdata()` produced `isotopeLabel.x`/`.y` suffixes because `isotopeLabel` was already parsed from row IDs. Fixed by excluding columns already present in `long_data` from the join.

### DESCRIPTION

Added `methods` to Imports (for `new("EList", ...)`). Added `@importFrom methods new` to `.lfqdata_to_elist`.

### Vignettes

- `LimmaBackend.Rmd` — added vooma section (`build_model_limma_voom`) and two limpa examples: (1) peptide → protein aggregation, (2) peptide-level impute-only analysis. Replaced `cat()` output with `knitr::kable()` tables throughout.
- `ContrastFacades.Rmd` — added limpa facade to the parallel comparison (protein-level facades, volcano plots, missing protein tables, rescued estimates).

### Tests

48 new tests in `tests/testthat/test-ContrastsLimpa.R`:
- `strategy_limpa` creation
- `AggregateLimpa` aggregation (SE/n_obs columns, no NAs, hierarchy reduction, wide pivotable)
- `AggregateLimpa` impute_only (hierarchy preserved, NAs filled)
- `build_model_limpa` returns `ModelLimma`
- `ContrastsLimma` works with limpa `ModelLimma`
- `ContrastsLimpaFacade` end-to-end (get_contrasts, get_missing, to_wide)
- `build_contrast_analysis(method = "limpa")`
- Rejection when `opt_se` is not set

### Files modified

| File | Changes |
|------|---------|
| `R/AnalysisConfiguration.R` | Added `opt_se` field (3 lines) |
| `R/LFQDataAggregator.R` | Added `AggregateLimpa`, `.elist_to_lfqdata()` (~180 lines) |
| `R/ContrastsLimpa.R` | **New file**: `StrategyLimpa`, `strategy_limpa()`, `build_model_limpa()` |
| `R/ContrastsLimma.R` | Extracted `.lfqdata_to_elist()`, `.resolve_weights()`; refactored 5 `build_model_*` functions |
| `R/ContrastsFacades.R` | Added `ContrastsLimpaFacade` |
| `R/build_contrast_analysis.R` | Added `"limpa"` method |
| `DESCRIPTION` | Added `methods` to Imports |
| `NAMESPACE` | Added limpa exports + `importFrom(methods, new)` |
| `tests/testthat/test-ContrastsLimpa.R` | **New file**: 48 tests |
| `tests/testthat/test-ContrastsLimma.R` | Added vooma/voom_impute test cases |
| `vignettes/LimmaBackend.Rmd` | Added vooma + two limpa examples, `cat()` → `knitr::kable()` |
| `vignettes/ContrastFacades.Rmd` | Added limpa facade to comparison |

All 618 tests pass (0 failures, 42 pre-existing warnings).

---

## 2026-03-26 — Remove Unused LFQDataImp R6 Class

Removed `LFQDataImp` R6 class and `LFQData$get_Imputer()` factory method. The class was a thin wrapper around `impute_with_zcomp()` — defined and tested but never used in any workflow, vignette, or downstream package (prolfquapp, prophosqua, prolfquabenchmark).

The standalone exported functions (`estimate_lod_global`, `function_lod_quantile`, `impute_with_zcomp`) remain in `R/LFQDataImp.R`. Documentation cross-references in `LFQDataPlotter` and `tidyMS_plotting.R` updated to point to `impute_with_zcomp()`.

**Files modified:** `R/LFQDataImp.R`, `R/LFQData.R`, `R/LFQDataPlotter.R`, `R/tidyMS_plotting.R`, `tests/testthat/test-LFQData.R`. All 548 tests pass.

---

## 2026-03-26 — Refactor LFQDataAggregator into Strategy Classes

Replaced monolithic `LFQDataAggregator` (one class, four methods) with three focused strategy classes sharing a uniform `$aggregate()` API:

| Class | Replaces | Input requirement |
|-------|----------|-------------------|
| `AggregateMedpolish` | `$medpolish()` | Transformed (log) intensities |
| `AggregateRlm` | `$lmrob()` | Transformed (log) intensities |
| `AggregateTopN` | `$sum_topN()` / `$mean_topN()` | Raw intensities; `func = "sum"` or `"mean"` |

**Design:** Composition, not inheritance. Shared logic (`plot`, `write_plots`) extracted to package-internal helpers (`.aggregator_plot()`, `.aggregator_write_plots()`). Validation warnings moved from method calls to constructors.

**Factory updated:** `LFQData$get_Aggregator(method = "medpolish", ...)` dispatches to the correct class. `...` passes `prefix`, `N`, `func` to the constructor.

**Bug fix:** `aggregate_intensity_topN()` used `match.call()` to capture `N` for column naming. When `N` was passed as an expression (e.g. `self$N` from an R6 method), `match.call()` captured the unevaluated expression, causing `:=` to fail with newer rlang. Fixed by using `N` directly instead of `xcall$N`.

### Files modified

| File | Changes |
|------|---------|
| `R/LFQDataAggregator.R` | Replaced `LFQDataAggregator` with `AggregateMedpolish`, `AggregateRlm`, `AggregateTopN` + shared helpers |
| `R/LFQData.R` | `get_Aggregator(method, ...)` factory with `switch` dispatch |
| `R/tidyMS_aggregation.R` | Removed fragile `match.call()` in `aggregate_intensity_topN()` |
| `R/tidyMS_reshaping.R` | Updated roxygen example |
| `tests/testthat/test-LFQData.R` | Split into 3 test blocks (one per strategy class) |
| `tests/testthat/test-ContrastsFacades.R` | Updated to `AggregateMedpolish$new()$aggregate()` |
| `vignettes/ContrastFacades.Rmd` | Updated aggregation calls |
| `vignettes/ContrastFacade2Factor.Rmd` | Updated aggregation calls |
| `vignettes/SimulateData.Rmd` | Updated aggregation calls |

### Downstream consumer updates (prolfquapp)

| File | Changes |
|------|---------|
| `prolfquapp/R/aggregation_IBAQ.R` | `get_Aggregator(method)` + `$aggregate()` |
| `prolfquapp/R/R6_ProteinDataPrep.R` | `get_Aggregator(method)` + `$aggregate()` |

All 549 tests pass (prolfqua). All 6 tests pass (prolfquapp).

---

## 2026-03-26 — Limma LOD Imputation Backend + lm_impute nr_children Fix

### `ContrastsLimmaImputeFacade` and `build_model_limma_impute()`

From `TODO/TODO_limma_with_LOD.md`. Added LOD imputation for limma, mirroring the existing `lm_impute` approach but adapted for limma's matrix-based fitting.

**New functions (`R/ContrastsLimma.R`):**
- `compute_borrowed_variance_limma()` — computes median sigma and df from successful proteins in an `MArrayLM` fit
- `build_model_limma_impute()` — fits limma on original data, identifies failed proteins (NA coefficients), imputes NAs with LOD, refits, splices imputed rows back into the fit with borrowed sigma and corrected df

**Key design decisions:**
- Only sigma borrowing (no `borrow_method`). Limma represents variance as `sigma^2 * stdev.unscaled^2`; replacing sigma alone is sufficient
- Direct field replacement on `MArrayLM` (coefficients, stdev.unscaled, sigma, Amean, df.residual) — no new S3 class needed
- df correction: imputed proteins use `max(n_observed - p, 1)` (default) or borrowed median df, never the inflated df from the LOD-imputed fit

**New facade (`R/ContrastsFacades.R`):**
- `ContrastsLimmaImputeFacade` — same API as all other facades (`get_contrasts()`, `get_missing()`, `get_Plotter()`, `to_wide()`)

**Also added:** `ContrastsLimma(eBayes = FALSE)` parameter to return raw unmoderated statistics, enabling downstream DEqMS moderation on limma results.

### `lm_impute` nr_children fix

`lm_impute` returned 157 results instead of the expected 160. Root cause: `.impute_one_protein()` joins with a sample template to add rows for missing groups, but the template only contains annotation columns — so `nr_children` was NA on new rows. Since `strategy_lm` passes `weights = nr_children` to `lm()`, the NA weights caused the refit to fail.

Fix: thread `nr_children_col` from `build_model_impute()` → `impute_refit_singular()` → `.impute_one_protein()`, and after the join, fill NA `nr_children` with 1 (conservative weight for imputed rows with no observed peptides).

### Vignettes

- `ContrastFacades.Rmd` — updated with `limma_impute` and `lm_impute` facades, split volcano plots (standard vs imputation with rescued protein highlighting), missing protein tables, rescued estimates comparison
- `ContrastFacade2Factor.Rmd` — restructured to match `ContrastFacades.Rmd`: added `lm_impute` and `limma_impute`, rescued protein flagging, split volcanos with red diamond markers, model name counts, missing protein tables, per-sample intensities, rescued estimates

### Files modified

| File | Changes |
|------|---------|
| `R/ContrastsLimma.R` | `eBayes` param on `ContrastsLimma`, `compute_borrowed_variance_limma()`, `build_model_limma_impute()` |
| `R/ContrastsFacades.R` | `ContrastsLimmaImputeFacade`, updated `FACADE_REGISTRY` |
| `R/build_contrast_analysis.R` | Added `"limma_impute"` method |
| `R/tidyMS_build_model.R` | `nr_children_col` parameter threaded through imputation chain; `annotation_vars()` sample template |
| `tests/testthat/test-ContrastsLimma.R` | 6 new test cases for limma impute |
| `vignettes/ContrastFacades.Rmd` | Full restructure with imputation facades and rescued protein visualization |
| `vignettes/ContrastFacade2Factor.Rmd` | Aligned structure with `ContrastFacades.Rmd` |

All 547 tests pass.

---

## 2026-03-25 — Performance Review Items 8, 9, 10

From `TODO/TODO_perforance_review.md`. Final three items — all 10 now complete.

- **Item 8 (`apply()` → `rowSums` + partial sort):** `estimate_lod_global()` replaced `apply(data_matrix, 1, function(x) sum(is.na(x)))` with vectorized `rowSums(is.na(data_matrix))`. `function_lod_quantile()` no longer converts matrix to data.frame; uses `sort(x, partial = n)[1:n]` for O(n) partial sort instead of full `sort()`. (`R/LFQDataImp.R`)

- **Item 9 (`get_LOD()` caching):** Added `private$.lod_cache` to `MissingHelpers` R6 class. `get_LOD()` now computes once on first call and returns cached value thereafter, matching the existing `self$stats` caching pattern. Eliminates redundant filter + quantile computation on every call. (`R/tidyMS_missingness_imputation.R`)

- **Item 10 (filter-compute-rejoin):** Replaced the filter → mutate → select → `left_join` pattern in `model_analyse()` with a single `mutate` using `purrr::map2_*` that checks `has_model_fit` per row. Failed fits get `NA` directly from the else branch, eliminating the intermediate filtered tibble and rejoin. (`R/tidyMS_build_model.R`)

All 534 tests pass.

---

## 2026-03-25 — Performance Review Item 7: Single pivot_wider in pivot_model_contrasts_2_Wide

From `TODO/TODO_perforance_review.md`.

The function looped over value columns, calling `pivot_wider` per column via the `m_spread` helper, then chained N-1 `left_join` calls via `purrr::reduce`. Replaced with a single `tidyr::pivot_wider(values_from = columns, names_glue = "{.value}.{contrast}")`. Eliminates the helper function, the loop, and the reduce/left_join chain.

**Files:** `R/tidyMS_contrasts.R`. All 534 tests pass.

---

## 2026-03-25 — Performance Review Items 5 & 6: Vectorized Contrast Computation

From `TODO/TODO_perforance_review.md`.

- **Item 5 (`compute_contrast` vectorization):** The original loops per-row of the `linfct` matrix, creating one `data.frame()` per contrast, then `bind_rows()`. The vectorized variant (`compute_contrast_vectorized`) uses matrix multiplication (`linfct %*% coef`) with zero-filled coefficients and invalid-row masking for NA coefficients. Standard errors computed via `sqrt(diag(L %*% Sigma %*% t(L)))` in one pass.

- **Item 6 (`linfct_matrix_contrasts` vectorization):** The original loops `rlang::parse_expr()` + `dplyr::mutate()` per contrast expression. The vectorized variant (`linfct_matrix_contrasts_vectorized`) pre-parses all expressions with `lapply(parse_expr)` and evaluates in a single `dplyr::mutate(data, !!!parsed)` call. Falls back to per-expression loop on error for granular reporting.

### Hot-swap mechanism

Both originals are **untouched**. Vectorized variants live alongside them in `R/tidyMS_contrasts.R`. A package option controls dispatch:

```r
options(prolfqua.vectorize = TRUE)   # use vectorized path
options(prolfqua.vectorize = FALSE)  # use original (default)
```

`compute_contrast()` and `linfct_matrix_contrasts()` check `getOption("prolfqua.vectorize")` at the top and delegate accordingly. This affects all Wald test facades (lm, rlm, firth, lmer) and limma's `linfct_matrix_contrasts` usage — no changes to R6 classes or strategy constructors needed.

### Test harness

`tests/testthat/test-vectorize-contrasts.R` (13 tests, 49 assertions):
- **Section A** (5 tests): `linfct_matrix_contrasts` vs `_vectorized` — parallel3, self-referencing factors, interaction model, single contrast, invalid contrast warning
- **Section B** (7 tests): `compute_contrast` vs `_vectorized` — full model, parallel3, single-row, all pairwise, different confint, manual NA coefficients, simulated incomplete models
- **Section C** (1 test): End-to-end facade hot-swap via `build_contrast_analysis(method = "lm")` with `options(prolfqua.vectorize)` toggle — verifies diff, std.error, p.value, FDR match to 1e-12

### Documentation

Roxygen `@section Vectorized mode` added to `build_contrast_analysis()`. One-liner notes added to `Contrasts`, `ContrastsLMFacade`, `ContrastsRLMFacade`, `ContrastsFirthFacade`, `ContrastsLimmaFacade`.

**Files:** `R/tidyMS_contrasts.R`, `R/build_contrast_analysis.R`, `R/Contrasts.R`, `R/ContrastsFacades.R`, `tests/testthat/test-vectorize-contrasts.R`. All 529 tests pass.

---

## 2026-03-25 — Performance Review Items 1–4

From `TODO/TODO_perforance_review.md`.

- **Item 1 (O(n²) copy-on-modify):** `impute_refit_singular` modified tibble list-columns element-by-element in a for-loop, triggering copy-on-modify per iteration. Fixed by collecting results into a list first, then bulk-assigning with `[idx] <-`. (`R/tidyMS_build_model.R`)

- **Item 2 (redundant `complete_cases`):** `summarize_stats_factors` called `summarize_stats` N+1 times, each calling `tidyr::complete()` (expensive cross-join) on the same data. Added `.completed` parameter to `summarize_stats`/`summarize_stats_all`; callers now complete once and pass `.completed = TRUE`. Also fixed the `"everything"` branch in `LFQDataStats`. (`R/tidyMS_stats.R`, `R/LFQDataStats.R`)

- **Item 3 (sequential mutates):** Combined 5 sequential `mutate()` calls into one in `model_analyse` and 3 sequential `mutate()` calls into one in `moderated_p_limma`. Each avoided mutate was copying a tibble with nested model objects. (`R/tidyMS_build_model.R`, `R/tidyMS_moderation.R`)

- **Item 4 (decorator deep-cloning):** Original claim was mostly wrong — critical review showed all 4 cloning decorators (Transformer, Aggregator, Stats, Imp) are justified because they all mutate data or config. `LFQDataSummariser` already correctly avoids cloning. **Bug found and fixed:** `LFQDataPlotter` did NOT clone but DID mutate the original data via `na.omit()`, silently dropping NA rows from the caller's `LFQData`. Fixed by adding `$clone(deep = TRUE)`. (`R/LFQDataPlotter.R`)

All 485 tests pass after each change.

---

## 2026-03-25 — Item 9: Split `AnalysisConfiguration.R` into cohesive files

From `TODO/TODO_10_refactorings.md` Item 9. Split the 870-line file (15 exports, 3 concerns) into cohesive files based on caller analysis:

| Destination | Functions moved | Primary callers |
|---|---|---|
| `R/tidyMS_data_setup.R` | `setup_analysis`, `complete_cases`, `sample_subset`, `separate_hierarchy` | LFQData.R, simulate_LFQ_data.R |
| `R/tidyMS_summarize_hierarchy.R` | `table_factors`, `table_factors_size`, `hierarchy_counts`, `HierarchyCountsSample`, `hierarchy_counts_sample`, `summarize_hierarchy` | LFQData.R, LFQDataSummariser.R, LFQDataStats.R |
| `R/utilities.R` | `make_interaction_column` | tidyMS_contrasts.R, tidyMS_build_model.R, tidyMS_plotting.R, tidyMS_stats.R |

`AnalysisConfiguration.R` reduced from 870 to ~345 lines (config container + serialization only). No API changes — all exports preserved, NAMESPACE unchanged. All 485 tests pass.

---

## 2026-03-25 — Item 10: Deduplicate `contrasts_linfct` branches

From `TODO/TODO_10_refactorings.md` Item 10. The function had two near-identical code paths for matrix vs list `linfct` inputs (only difference: `linfct` vs `linfct[[i]]`). Normalized the input at the top — if matrix, replicate as a list — then a single loop handles both cases. Net -7 lines.

**Files:** `R/tidyMS_contrasts.R`. All 485 tests pass.

---

## 2026-03-25 — Items 6 & 7: Extract helpers from `tidyMS_build_model.R`

From `TODO/TODO_10_refactorings.md` Items 6 and 7.

- **Item 6:** Extracted `.impute_one_protein()` from `impute_refit_singular` for-loop body. The helper takes one protein's data, imputes missing values with LOD, refits the model, and returns a named list (or `NULL` on failure). The orchestration loop remains in `impute_refit_singular`.
- **Item 7:** Extracted `plot_lrt_diagnostics()` from `LR_test`. The plotting code (histogram + ECDF PDF) is now a separate internal function. `LR_test` calls it when `path` is non-NULL.

Also renamed PascalCase local variables to snake_case across Model classes: `Model_Coeff` → `model_coeff`, `Model_Anova` → `model_anova`, `VolcanoPlot` → `volcano_plot` in `Model.R`, `ModelFirth.R`, `ContrastsLimma.R`.

**Files:** `R/tidyMS_build_model.R`, `R/Model.R`, `R/ModelFirth.R`, `R/ContrastsLimma.R`. All 485 tests pass.

---

## 2026-03-25 — Item 5: Rename `exists_lmer` column → `has_model_fit`

From `TODO/TODO_10_refactorings.md` Item 5. Renamed the `exists_lmer` column to `has_model_fit` across all R source, tests, and vignettes. The column is a boolean indicating whether a model fit succeeded (not lmer-specific), so the new name is accurate for all backends (lm, rlm, lmer, logistf).

**Files:** `R/tidyMS_build_model.R`, `R/Contrasts.R`, `R/tidyMS_R6_Modelling.R`, `R/simulate_LFQ_data.R`, `tests/testthat/test-ImputeModel.R`, `vignettes/TestingMissingInference.Rmd`. All 485 tests pass.

---

## 2026-03-25 — Item 4: Rename contrast functions and consolidate `.ehandler`

From `TODO/TODO_10_refactorings.md` Item 4. Renamed exported contrast functions and consolidated duplicated error handler:

- `my_contrast` → `.compute_contrast` (made non-exported; only used internally by `compute_contrast` and `contrasts_linfct_firth`)
- `my_contrast_V2` → `compute_contrast` (sole exported lm/rlm/logistf contrast function)
- `my_contest` → `compute_lmer_contrast` (lmer-specific contrast wrapper)
- `.ehandler` → `.error_handler` (consolidated from 2 identical copies in `tidyMS_R6_Modelling.R` and `tidyMS_plotting.R` into single definition in `utilities.R`)

**Files:** `R/tidyMS_contrasts.R`, `R/tidyMS_R6_Modelling.R`, `R/logistf.R`, `R/tidyMS_plotting.R`, `R/utilities.R`, `R/tidyMS_moderation.R`. All 485 tests pass.

---

## 2026-03-25 — Item 3: Standardize `nrcoeff` → `nr_coef` / `nr_coef_not_NA`

From `TODO/TODO_10_refactorings.md` Item 3. Standardized column names to single-f snake_case (`nr_coef`, `nr_coef_not_NA`) to match R's own `coef()` convention. Renamed closures: `count_coefficients` → `count_coef`, `nrcoeff_not_NA` → `count_coef_not_NA`. Renamed local var `mcoef` → `max_coef` in Contrasts.R and ContrastFirth.R.

**Files:** `R/tidyMS_build_model.R`, `R/Model.R`, `R/ModelFirth.R`, `R/Contrasts.R`, `R/ContrastFirth.R`, `vignettes/TestingMissingInference.Rmd`. All 485 tests pass, all vignettes build.

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
