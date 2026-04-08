# Plan: Phases 2–4 — function parameters, function names, R6 methods

## Context

Phase 1 (local variable renames) committed. Now three phases, all hard-switch, prolfqua + downstream. Data frame column name strings stay as-is.

---

## Phase 2: Function parameter renames

All calls are positional in downstream — no downstream breakage. Rename in definitions + fix internal named call sites.

| Old param | New param | Files (definition) |
|-----------|-----------|-------------------|
| `colName` | `col_name` | AnalysisConfiguration.R:78 `set_response()` |
| `workIntensity` | `work_intensity` | AnalysisConfiguration.R:266 `make_reduced_hierarchy_config()` |
| `modelName` | `model_name` | ContrastsLimma.R:228,913; ContrastsLimpa.R:84; ContrastsInterface.R:76; ContrastsROPECA.R:66; tidyMS_build_model.R:84 |
| `modelName_Int` | `model_name_int` | tidyMS_build_model.R:84 `plot_lrt_diagnostics()` |
| `preserveMean` | `preserve_mean` | LFQDataTransformer.R:123,134,241,377,464 |
| `sampleName` | `sample_name` | tidyMS_aggregation.R:233,291 |
| `factorDepth` | `factor_depth` | LFQData.R:145 `omit_NA()` |
| `nrNA` | `nr_na` | LFQData.R:145 `omit_NA()` |
| `proteinName` | `protein_name` | tidyMS_aggregation.R:110 `plot_hierarchies_line()` |
| `intesityNewName` | `intensity_new_name` | LFQDataTransformer.R:241 (fix typo too) |
| `modelDF` | `model_df` | tidyMS_build_model.R:365 `compute_borrowed_variance()` |
| `modelProteinF` | `complete_models` | tidyMS_build_model.R:584 `get_complete_model_fit()` |
| `modelProteinF` / `modelProteinF_Int` | `complete_models` / `complete_models_int` | tidyMS_build_model.R:39-49 `compare_models_lrt()` |
| `modelWithInteractionsContrasts` | `model_interaction_contrasts` | tidyMS_contrasts.R:626 `pivot_model_contrasts_2_Wide()` |

### Internal named call sites to fix

| Call site | Change |
|-----------|--------|
| LFQDataTransformer.R:124 | `preserveMean = preserveMean` → `preserve_mean = preserve_mean` |
| LFQDataTransformer.R:140 | same |
| tidyMS_aggregation.R:292-293,327-331 | `sampleName = sampleName` → `sample_name = sample_name` |
| tidyMS_aggregation.R:120 | `proteinName = proteinName` → `protein_name = protein_name` |

### Inner function renames (local helper functions)

| Old | New | File |
|-----|-----|------|
| `getCoeffs` | `get_coeffs` | tidyMS_contrasts.R:47 |
| `getSampleSize` | `get_sample_size` | tidyMS_stats.R:326,425,470 |
| `getAST` | `get_ast` | tidyMS_missingness_summary.R:20 |

### Roxygen `@param` tags

Update all `@param` roxygen tags that reference the old parameter names.

**Commit after Phase 2. Verify: `make check-fast`.**

---

## Phase 3: Exported function name renames

| Old | New | File | Callers (prolfqua) | Callers (downstream) |
|-----|-----|------|--------------------|---------------------|
| `scriptCopyHelperVec` | `script_copy_helper_vec` | vignetteHelpers.R:23 | — | prolfquapp/R/copy_helpers.R (3×), prolfquasaint/R/copy_helpers.R (2×), prophosqua/R/prophosqua_copy_helpers.R (1×) |
| `pivot_model_contrasts_2_Wide` | `pivot_model_contrasts_to_wide` | tidyMS_contrasts.R:626 | ContrastsLimma.R:1078, Contrasts.R:254 | prolfquasaint/R/ContrastSaintExpress.R:116 |
| `isSingular_lm` | `is_singular_lm` | tidyMS_build_model.R:564 | tidyMS_R6_Modelling.R:200 | — |
| `sim_lfq_data_2Factor_config` | `sim_lfq_data_2factor_config` | simulate_LFQ_data.R:268 | logistf.R, tidyMS_stats.R, tidyMS_missingness_imputation.R, tidyMS_build_model.R, tests (3), vignettes (3) | — |
| `UpSet_interaction_missing_stats` | `upset_interaction_missing_stats` | tidyMS_missingness_summary.R:247 | LFQDataSummariser.R | prolfquapp/R/R6_DEAReportGenerator.R, prolfquasaint/inst/ (2 files) |
| `get_UniprotID_from_fasta_header` | `get_uniprot_id_from_fasta_header` | utilities.R:23 | — | prolfquapp/R/R6_ProteinAnnotation.R:485, prolfquasaint/inst/ (3×) |

Also update NAMESPACE (via `make document`) and any roxygen `@examples` that call these functions.

**Commit after Phase 3. Verify: `make check-fast` in prolfqua, then grep downstream packages and fix callers, test downstream.**

---

## Phase 4: R6 method renames

### 4a: Delete deprecated aliases (fix callers first)

| Alias | Replacement | Callers to fix |
|-------|------------|----------------|
| `hierarchyKeys()` | `hierarchy_keys()` | prolfquapp/R/R6_DEAReportGenerator.R:457, prolfquapp/R/R6_QC_Abundances.R:49 |
| `hkeysDepth()` | `hierarchy_keys_depth()` | no current callers |

Then delete both from AnalysisConfiguration.R (lines ~122-124, ~143-145).

### 4b: Rename R6 methods

| Old | New | Class | File | Callers (prolfqua) | Callers (downstream) |
|-----|-----|-------|------|--------------------|---------------------|
| `subject_Id()` | `subject_id()` | LFQData | LFQData.R:106 | ContrastsFacades.R (4×), logistf.R (2×), tidyMS_build_model.R (2×) | — |
| `omit_NA()` | `omit_na()` | LFQData | LFQData.R:145 | examples only (LFQData.R:31,34) | — |
| `NA_heatmap()` | `na_heatmap()` | LFQDataPlotter | LFQDataPlotter.R:165 | — | prolfquapp vignettes (2×), prolfquasaint/inst/ (1×) |
| `get_LOD()` | `get_lod()` | MissingHelpers | tidyMS_missingness_imputation.R:68 | tidyMS_missingness_imputation.R (2×), tidyMS_build_model.R:203, ContrastsLimma.R:322,582 | — |

**Commit after Phase 4. Verify: `make check-fast` in prolfqua, fix + test downstream.**

---

## Execution order

1. Phase 2 — parameter renames (prolfqua only) → commit
2. Phase 3 — function name renames (prolfqua + downstream) → commit
3. Phase 4a — delete deprecated aliases (prolfqua + prolfquapp) → commit
4. Phase 4b — R6 method renames (prolfqua + downstream) → commit

Each commit verified with `make check-fast`.
