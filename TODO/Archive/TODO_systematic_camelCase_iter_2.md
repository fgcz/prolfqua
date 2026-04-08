# Plan: snake_case migration + R6 encapsulation for prolfqua

## Context

The `AnalysisConfiguration` field migration is done. Now we address camelCase **throughout** prolfqua and simultaneously improve R6 encapsulation: minimize direct field access, expose methods instead. Hard-switch strategy (no deprecation wrappers). Data frame column names stay as-is.

## Scope

| In scope | Out of scope |
|----------|-------------|
| All variable names (local, parameters, fields, methods) in prolfqua | Data frame column name strings (`"sampleName"`, `"isotopeLabel"`) |
| R6 field encapsulation (private + accessors) | prolfquabenchmark (deferred) |
| Downstream callers in prolfquapp, prophosqua, prolfquasaint, prolfquappPTMreaders | Factory method names `get_Transformer()` etc. (these reference class names — keep) |

---

## Phase 1: Local variables inside function bodies (prolfqua only)

**Why first:** Zero API impact, pure cleanup, builds momentum. Do file-by-file.

### Rename table for common local variables

| Old | New | Files affected |
|-----|-----|----------------|
| `valueVars` | `value_vars` | AnalysisConfiguration.R, LFQData.R |
| `annotationVars` | `annotation_vars` | AnalysisConfiguration.R |
| `idVars` / `id_vars` | already snake_case in some places | tidyMS_data_setup.R |
| `colName` (local, not param) | `col_name` | various |
| `factorDepths` | `factor_depths` | tidyMS_contrasts.R |
| `pdata` | keep — it's a convention for "prepared data" | — |

### Files to sweep (alphabetical)

Each file: read, rename camelCase locals, verify no column-name strings were changed.

- `R/AnalysisConfiguration.R` — `valueVars`, `annotationVars` in `value_vars()` / `annotation_vars()` methods
- `R/Contrasts.R` — locals in methods
- `R/ContrastsFacades.R` — locals in facade methods
- `R/ContrastsLimma.R` — locals
- `R/ContrastsPlotter.R` — locals
- `R/LFQData.R` — locals in `nr_B_in_A()`, `filter_proteins_by_peptide_count()`
- `R/LFQDataAggregator.R` — `workIntensity`, `hierarchyDepth` as local strings
- `R/LFQDataTransformer.R` — locals in `transform_work_intensity()`
- `R/Model.R` — locals
- `R/tidyMS_aggregation.R` — `sampleName` locals (careful: some are column name strings)
- `R/tidyMS_build_model.R` — locals
- `R/tidyMS_contrasts.R` — `factorDepths`, `modelWithInteractionsContrasts`
- `R/tidyMS_moderation.R` — locals
- `R/tidyMS_stats.R` — locals
- `R/tidyMS_plotting.R` — locals
- `R/utilities.R` — locals

**Verification:** `make test` after each batch of ~3-4 files.

---

## Phase 2: Exported function parameter names (prolfqua + downstream)

**Hard switch.** Rename parameters, then grep-fix all callers in downstream packages.

| Parameter | Functions using it | New name |
|-----------|-------------------|----------|
| `colName` | `set_response()` | `col_name` |
| `workIntensity` | `make_reduced_hierarchy_config()` | `work_intensity` |
| `modelName` | `build_model_limma()`, `build_model_limpa()`, `merge_contrasts_results()`, `plot_lrt_diagnostics()`, ContrastsLimma init, ContrastsROPECA init | `model_name` |
| `preserveMean` | `robust_scale()`, `robscale()`, `robscale_subset()`, `scale_with_subset()`, `transform_work_intensity()` | `preserve_mean` |
| `sampleName` | `medpolish_estimate()`, `medpolish_estimate_df()` | `sample_name` |
| `factorDepth` | `omit_NA()` | `factor_depth` |
| `proteinName` | `plot_hierarchies_line()` | `protein_name` |
| `intesityNewName` | `transform_work_intensity()` | `intensity_new_name` (fix typo too) |
| `modelDF` | `compute_borrowed_variance()` | `model_df` |
| `modelName_Int` | `plot_lrt_diagnostics()` | `model_name_int` |
| `modelWithInteractionsContrasts` | `pivot_model_contrasts_2_Wide()` | `model_interaction_contrasts` |

**Downstream callers to fix:** Search each parameter name in prolfquapp, prophosqua, prolfquasaint, prolfquappPTMreaders.

**Verification:** `make test` in prolfqua, then `make test` in each downstream package.

---

## Phase 3: Exported function names (prolfqua + downstream)

Hard switch. Rename, grep-fix all callers.

| Old | New | File |
|-----|-----|------|
| `scriptCopyHelperVec` | `script_copy_helper_vec` | vignetteHelpers.R |
| `pivot_model_contrasts_2_Wide` | `pivot_model_contrasts_to_wide` | tidyMS_contrasts.R |
| `get_UniprotID_from_fasta_header` | Move to prolfquapp or rename `get_uniprot_id_from_fasta_header` | utilities.R |
| `isSingular_lm` | `is_singular_lm` | tidyMS_build_model.R |
| `sim_lfq_data_2Factor_config` | `sim_lfq_data_2factor_config` | simulate_LFQ_data.R |
| `UpSet_interaction_missing_stats` | `upset_interaction_missing_stats` | tidyMS_summarize_hierarchy.R |

**Decision needed:** `get_UniprotID_from_fasta_header` — user wants to drop or move to prolfquapp. Check callers first.

---

## Phase 4: R6 method names (prolfqua + downstream)

### 4a: Delete deprecated aliases

| Method | Replacement | Callers to fix |
|--------|------------|----------------|
| `hierarchyKeys()` | `hierarchy_keys()` | `LFQData.R:385`, prolfquapp `R6_DEAReportGenerator.R:457`, prolfquapp `R6_QC_Abundances.R:49` |
| `hkeysDepth()` | `hierarchy_keys_depth()` | `LFQData.R:384` |

Fix callers, then delete from AnalysisConfiguration.R.

### 4b: Rename other camelCase methods

| Method | Class | New name | Callers (approx) |
|--------|-------|----------|-----------------|
| `subject_Id()` | `LFQData` | `subject_id()` | internal (~5), prolfquapp (~3) |
| `omit_NA()` | `LFQData` | `omit_na()` | internal (~2), prolfquapp (~3) |
| `NA_heatmap()` | `LFQDataPlotter` | `na_heatmap()` | internal (~1), prolfquapp (~2) |
| `get_LOD()` | `MissingHelpers` | `get_lod()` | internal only (~2) |

### 4c: DO NOT rename (class name references)

These stay because the capitalized part is a **class name**:
- `get_Transformer()`, `get_Aggregator()`, `get_Plotter()`, `get_Stats()`, `get_Summariser()`, `get_Imputer()`

---

## Phase 5: R6 field names — rename to snake_case

| Old field | New field | Classes | Downstream callers |
|-----------|-----------|---------|-------------------|
| `modelName` | `model_name` | Contrasts, ContrastsPlotter, ContrastsTable, ContrastsMissing, ContrastsModerated, ContrastsFirth, ContrastsROPECA, Model, ModelFirth | prolfquasaint `ContrastSaintExpress.R` (2) |
| `subject_Id` | `subject_id` | same set + ContrastsPlotter | prolfquasaint `ContrastSaintExpress.R` (2) |
| `contrastDF` | `contrast_df` | ContrastsPlotter | internal only |
| `modelDF` | `model_df` | Model | internal only |

**CRITICAL: `modelName` duality.** `modelName` is both an R6 field AND a data frame column name in output. The field rename (`$modelName` → `$model_name`) does NOT change the column name in output data frames. Verify all `get_contrasts()` / `get_coefficients()` methods emit column named `modelName` in their output — that string stays.

**Leave as-is:** `p.adjust`, `avg.abundance` (R stats convention with `.` separator).

---

## Phase 6: R6 encapsulation — make fields private, add accessors

This is the architectural improvement. Goal: external code accesses R6 objects only through methods, never through `$field` directly. This enables future internal representation changes.

### 6a: AnalysisConfiguration — add setter methods for structured fields

Currently `config$hierarchy[["protein_Id"]] <- "col"` is used ~25 times in downstream packages. Add methods:

```r
# New methods to add:
add_hierarchy = function(name, columns) {
  self$hierarchy[[name]] <- columns
  invisible(self)
},
add_factor = function(name, column) {
  self$factors[[name]] <- column
  invisible(self)
}
```

Existing methods already cover: `set_response()`, `get_response()`, `pop_response()`, `hierarchy_keys()`, `hierarchy_keys_depth()`, `factor_keys()`, `factor_keys_depth()`.

Still need getters/setters for: `file_name`, `sample_name`, `hierarchy_depth`, `factor_depth`, `isotope_label`, `ident_q_value`, `ident_score`, `is_response_transformed`, `nr_children`, `norm_value`.

**Decision point:** Do we add individual getters/setters for every scalar field, or keep scalar fields public and only encapsulate structured fields (`hierarchy`, `factors`, `work_intensity` stack)? Individual getters/setters for 15+ scalar fields is a lot of boilerplate.

**Recommendation:** Keep scalar config fields public for now. Encapsulate only:
- `hierarchy` (structured, list with named entries) → `add_hierarchy()`, `get_hierarchy()`
- `factors` (structured, list) → `add_factor()`, `get_factors()`
- `work_intensity` (stack semantics) → already has `set_response()` / `get_response()` / `pop_response()`

### 6b: LFQData — encapsulate `$data` and `$config`

Currently `lfqdata$data` is read 105+ times and written ~7 times in downstream packages. `lfqdata$config` is read 100+ times.

**New accessors:**
```r
# Already exists as methods — just need to formalize:
get_data = function() self$data,
get_config = function() self$config,
```

**Problem:** Downstream code does `lfqdata$data <- lfqdata$data |> dplyr::filter(...)` — direct mutation of the data field. Making `$data` private requires a `set_data()` or `update_data()` method, or the filter operations move into LFQData methods.

**Phased approach:**
1. First: Add `get_data()` / `get_config()` methods (non-breaking)
2. Migrate downstream callers from `lfqdata$data` → `lfqdata$get_data()`
3. Then: make `data` and `config` private

**Scale of downstream migration:**
- prolfquapp: ~205 accesses (`$data` + `$config`)
- prolfquasaint: ~29 accesses
- prophosqua: ~8 accesses
- prolfquappPTMreaders: ~19 accesses

### 6c: Decorator classes — `$lfq` field

The `$lfq` field on LFQDataTransformer etc. is the result access pattern:
```r
lfqdata <- lfqdata$get_Transformer()$log2()$robscale()$lfq
```

This is used ~1800 times across the ecosystem (including internal). Options:
1. **Keep `$lfq` public** — it's the standard decorator result access. Making it a method (`$get_lfq()`) adds verbosity for no real encapsulation benefit since the field is read-only in practice.
2. **Rename to `$result`** — more descriptive but massive blast radius.

**Recommendation:** Keep `$lfq` as-is. It's a well-established pattern and already snake_case.

### 6d: Contrasts/Model fields — encapsulate

| Field | Accessor to add | Notes |
|-------|----------------|-------|
| `model_name` (was `modelName`) | `get_model_name()` | Read-only after init |
| `subject_id` (was `subject_Id`) | `get_subject_id()` | Read-only after init |
| `contrast_df` (was `contrastDF`) | Already have `get_contrasts()` | `contrastDF` is internal storage, `get_contrasts()` is the public API |
| `model_df` (was `modelDF`) | Already have `get_coefficients()` / `get_anova()` | `modelDF` is internal storage |

**These fields can go private immediately** — external code should use `get_contrasts()`, `get_coefficients()`, etc.

---

## Phase 7: Refactor `function(data, config)` → `function(lfqdata)` signatures

The user specifically called this out. ~40 functions currently take separate `data` and `config` parameters. Many are only called from LFQData class methods that pass `self$data, self$config`.

### Category A: Functions called ONLY from LFQData methods (make private or internalize)

| Function | Called from | Action |
|----------|-----------|--------|
| `remove_small_intensities(pdata, config, ...)` | `LFQData$remove_small_intensities()` | Inline into class method |
| `filter_proteins_by_peptide_count(pdata, config)` | `LFQData$filter_proteins_by_peptide_count()` | Inline into class method |
| `complete_cases(pdata, config)` | `LFQData$complete_cases()` | Inline into class method |
| `nr_B_in_A(pdata, config)` | `LFQData$filter_difference()` | Inline or keep as internal helper |

### Category B: Functions called from BOTH class methods AND standalone (keep exported, change signature)

| Function | Current sig | New sig | Callers |
|----------|-----------|---------|---------|
| `summarize_stats(pdata, config, ...)` | `(pdata, config, ...)` | `(lfqdata, ...)` | LFQDataStats methods, standalone |
| `hierarchy_counts(pdata, config)` | `(pdata, config)` | `(lfqdata)` | LFQDataSummariser, standalone |
| `table_factors(pdata, config)` | `(pdata, config)` | `(lfqdata)` | LFQData$factors(), standalone |
| `plot_heatmap(data, config, ...)` | `(data, config, ...)` | `(lfqdata, ...)` | LFQDataPlotter, standalone |
| `plot_pca(data, config, ...)` | `(data, config, ...)` | `(lfqdata, ...)` | LFQDataPlotter, standalone |
| `transform_work_intensity(pdata, config, ...)` | `(pdata, config, ...)` | `(lfqdata, ...)` | LFQDataTransformer, standalone |
| `estimate_intensity(data, config, .func)` | `(data, config, .func)` | `(lfqdata, .func)` | LFQDataAggregator, standalone |

These functions would extract `data` and `config` from the lfqdata object internally:
```r
summarize_stats <- function(lfqdata, ...) {
  pdata <- lfqdata$data
  config <- lfqdata$config
  # ... rest unchanged
}
```

### Category C: Functions called standalone only (lower priority)

These are plotting/utility functions that take `(data, config)` but aren't wrapped by class methods. Change when touching the file for other reasons.

**Decision needed:** Is this phase worth the effort now, or defer until after Phases 1-6? It's a large refactor with moderate benefit — the main gain is API cleanliness.

---

## Phase 8: Fix remaining cross-package issues

### `AnalysisTableAnnotation` references
| File | Change |
|------|--------|
| `prolfquasaint/inst/application/SE2_DIANN/DIANN_SE.R:54` | `AnalysisTableAnnotation$new()` → `AnalysisConfiguration$new()` |
| `prolfquasaint/inst/application/SE2/CreateSaintExpress_Report.R:63` | Same |

### Stray camelCase config field access
| File | Change |
|------|--------|
| `prolfquapp/tests/testthat/Grp2Analysis_V2_R6.Rmd:114-115` | `ld$config$sampleName` → `ld$config$sample_name` |

---

## Execution Order

| Order | Phase | Scope | Breaking? | Effort |
|-------|-------|-------|-----------|--------|
| 1 | Phase 1 | Local variables | No | Medium (many files, zero risk) |
| 2 | Phase 4a | Delete `hierarchyKeys`/`hkeysDepth` aliases | Yes (4 callers) | Small |
| 3 | Phase 8 | Cross-package fixes | Yes (3 files) | Small |
| 4 | Phase 2 | Function parameters | Yes (downstream) | Medium |
| 5 | Phase 3 | Function names | Yes (downstream) | Small |
| 6 | Phase 4b | R6 method names | Yes (downstream) | Small |
| 7 | Phase 5 | R6 field names | Yes (downstream) | Medium |
| 8 | Phase 6 | Encapsulation (add accessors) | No (additive) | Medium |
| 9 | Phase 7 | `(data, config)` → `(lfqdata)` | Yes (downstream) | Large |

**Commit strategy:** Commit after each phase. Each phase should leave all tests passing.

## Verification per phase

```bash
cd /Users/wolski/projects/prolfqua_fml/prolfqua && make test
cd /Users/wolski/projects/prolfqua_fml/prolfquapp && make test
# After phases that touch downstream:
cd /Users/wolski/projects/prolfqua_fml/prophosqua && make test
cd /Users/wolski/projects/prolfqua_fml/prolfquasaint && make test
cd /Users/wolski/projects/prolfqua_fml/prolfquappPTMreaders && make test
```

## Open decisions

1. **`get_UniprotID_from_fasta_header`** — drop from prolfqua and move to prolfquapp? Need to check all callers first.
2. **Config scalar fields** — add individual getters/setters or keep public? Recommendation: keep public for now, encapsulate only structured fields.
3. **Phase 7 timing** — refactor `(data, config)` → `(lfqdata)` now or defer? It's the largest phase with moderate benefit.
4. **`$lfq` decorator field** — keep public (recommendation: yes).
