# Dead Code Analysis Report — prolfqua

## Context

Systematic dead code analysis across the prolfqua ecosystem (`prolfqua_fml/`). Verified against all downstream packages: prolfquapp, prophosqua, prolfquabenchmark, prolfquasaint, prolfquappPTMreaders. Out of ~204 exports, only ~32 functions/classes are directly referenced downstream.

---

## 1. `inst/deprecated/DEPRECATED.R` — ALL DEAD (delete entire file)

| Function | References | In NAMESPACE | Verdict |
|----------|-----------|-------------|---------|
| `aggregate_contrast` | None | No | Dead |
| `get_impute_contrasts_V1` | None | No | Dead |
| `get_imputed_contrasts` | None | No | Dead |
| `missigness_impute_factors_interactions` | None | No | Dead |
| `.missigness_impute_interactions` | None | No | Dead |

File is under `inst/`, not `R/` — never loaded. All contain `message("deprecated")`. Superseded by `ContrastsMissing` R6 class.

---

## 2. TIER 1: Definitely Unused — Safe to Remove (16 items)

These exported functions have **zero callers** anywhere in the ecosystem (no tests, no vignettes, no downstream packages).

| # | Function | File | Notes |
|---|----------|------|-------|
| 1 | `render_MQSummary_rmd` | `R/vignetteHelpers.R:112` | Vestigial MQ workflow |
| 2 | `spread_response_by_IsotopeLabel` | `R/AnalysisConfiguration.R:863` | Abandoned IsotopeLabel feature |
| 3 | `panel_hist` | `R/utilities.R:409` | Only in commented example |
| 4 | `pairs_w_abline` | `R/utilities.R:382` | Zero callers |
| 5 | `plot_heatmap_cor_iheatmap` | `R/tidyMS_plotting.R:329` | Uses iheatmapr, zero callers |
| 6 | `plot_screeplot` | `R/tidyMS_plotting.R:668` | Zero callers |
| 7 | `normalize_log2_robscale` | `R/tidyMS_R6_TransitionCorrelations.R:529` | Convenience wrapper, zero callers |
| 8 | `medpolish_protein_estimates` | `R/tidyMS_aggregation.R:960` | Dead (+ dead deps below) |
| 9 | `intensity_summary_by_hkeys` | `R/tidyMS_aggregation.R:847` | Only called by dead #8 |
| 10 | `response_as_matrix` | `R/tidyMS_aggregation.R:304` | Only called by dead #9 |
| 11 | `cor_jackknife_matrix` | `R/utilities.R:351` | Dead jackknife cluster |
| 12 | `cor_order` | `R/utilities.R:328` | Dead jackknife cluster |
| 13 | `jackknife` | `R/utilities.R:260` | Dead jackknife cluster |
| 14 | `jackknife_matrix` | `R/utilities.R:293` | Dead jackknife cluster |
| 15 | `scale_with_subset_by_factors` | `R/tidyMS_R6_TransitionCorrelations.R:485` | Zero callers |
| 16 | `nr_B_in_A_per_sample` | `R/tidyMS_R6_TransitionCorrelations.R:627` | Zero callers |
| 17 | `rank_by_NA` | `R/tidyMS_R6_TransitionCorrelations.R:821` | Zero callers |
| 18 | `create_config_MSFragger_MSstats` | `R/tidyMS_R6_ConcreteConfigurations.R` | Zero callers |
| 19 | `separate_factors` | `R/AnalysisConfiguration.R` | Zero callers |
| 20 | `plot_estimate` | `R/tidyMS_aggregation.R:698` | Zero callers |
| 21 | `my_contrast_V1` | `R/tidyMS_R6_Modelling.R:1007` | Superseded by V2 |

---

## 3. TIER 2: Probably Unused — Remove After Verification (18 items)

Used only internally or in roxygen examples. Not called by any downstream package.

| # | Function | File | Notes |
|---|----------|------|-------|
| 1 | `plot_hierarchies_line` | `R/tidyMS_aggregation.R` | Old aggregation plotting |
| 2 | `plot_hierarchies_line_df` | `R/tidyMS_aggregation.R` | Old aggregation plotting |
| 3 | `plot_raster` | `R/tidyMS_plotting.R:468` | Only via LFQDataPlotter |
| 4 | `plot_NA_heatmap` | `R/tidyMS_plotting.R:529` | Only via LFQDataPlotter |
| 5 | `filter_difference` | `R/tidyMS_MQ_workflow.R:53` | Only in one test |
| 6 | `create_config_Spectronaut_Peptide` | `R/tidyMS_R6_ConcreteConfigurations.R` | Only in data-raw script |
| 7 | `sample_subset` | `R/AnalysisConfiguration.R` | Only internal |
| 8 | `lfq_power_t_test_quantiles_V2` | `R/tidyMS_stats.R` | V2 suffix, one vignette |
| 9 | `LR_test` | `R/tidyMS_R6Model.R` | Only in one vignette |
| 10 | `matrix_to_tibble` | `R/utilities.R:222` | Only prolfquabenchmark |
| 11 | `interaction_missing_stats` | `R/tidyMS_missigness.R:47` | Deprecated, but still called by `UpSet_interaction_missing_stats` etc. |

---

## 4. De-export Candidates (keep function, remove `@export`)

Internal plumbing that works but should not be public API:

| Function | Reason |
|----------|--------|
| `isSingular_lm` | Internal strategy helper |
| `get_p_values_pbeta` | Internal |
| `contrasts_fisher_exact` | Internal |
| `model_analyse` | Internal build_model plumbing |
| `model_summary` | Internal, one vignette |
| `compute_pooled` | Called only from `poolvar` |
| `impute_with_zcomp` | Internal LFQDataImp |
| `function_lod_quantile` | Internal LFQDataImp |
| `estimate_lod_global` | Internal LFQDataImp |
| `which_missing` | Internal simulation helper |
| `moderated_p_deqms_long` | Internal ContrastsModeratedDEqMS |
| `panel_cor` | Internal helper for `pairs_smooth` |

---

## 5. TIER 3: Needs Your Decision

### Firth/GLM subsystem (`R/logistf.R`)
- `strategy_logistf`, `build_model_logistf`, `build_model_glm_peptide`, `build_model_glm_protein`, `ContrastsFirth`, `ContrastsFirthFacade`, `ModelFirth`, `contrasts_linfct_firth`
- Only used in prolfqua tests. Zero downstream usage. **Is anyone using Firth logistic regression?**

### ROPECA subsystem
- `ContrastsROPECA`, `ContrastsROPECAFacade`
- Only in prolfquabenchmark vignettes. **Keep for paper reproducibility?**

### Facade classes not referenced downstream
- `ContrastsDEqMSFacade`, `ContrastsLMFacade`, `ContrastsLMImputeFacade`, `ContrastsLmerFacade`, `ContrastsRLMFacade`
- **Need to verify against `FACADE_REGISTRY` contents**

---

## 6. Redundancy / Refactoring Opportunities

### `tidyMS_missigness.R` vs `tidyMS_missigness_V2.R`
- `interaction_missing_stats` is deprecated in favor of `summarize_stats_factors`
- But `missigness_histogram`, `missingness_per_condition`, `UpSet_interaction_missing_stats` still call the deprecated function
- **Should migrate callers to `summarize_stats_factors`**

### Dead aggregation cluster in `tidyMS_aggregation.R`
- `medpolish_protein_estimates` → `intensity_summary_by_hkeys` → `response_as_matrix` — entire chain is dead
- `plot_estimate`, `plot_hierarchies_line`, `plot_hierarchies_line_df` — old viz replaced by R6 classes

---

## 7. Downstream Usage Baseline (for reference)

### Actively used R6 classes (9)
`AnalysisConfiguration`, `LFQData`, `Contrasts`, `ContrastsModerated`, `ContrastsMissing`, `ContrastsInterface`, `ContrastsPlotter`, `ContrastsLMMissingFacade`, `ContrastsLimmaFacade`

### Actively used functions (~32)
`setup_analysis`, `find_package_file`, `R6_extract_values`, `scriptCopyHelperVec`, `nr_obs_experiment`, `nr_obs_sample`, `get_UniprotID_from_fasta_header`, `annotation_add_contrasts`, `UpSet_interaction_missing_stats`, `FACADE_REGISTRY`, `volcano_plotly`, `scatter_plotly`, `table_facade`, `sim_lfq_data_peptide_config`, `sim_lfq_data_protein_config`, `sim_lfq_data_2Factor_config`, `old2new`, `prolfqua_data`, `strategy_lm`, `build_model`, `merge_contrasts_results`, `list_to_AnalysisConfiguration`, `create_config_MQ_peptide`, `LFQDataToSummarizedExperiment`, `matrix_to_tibble`, `adjust_p_values`, `summary_ROPECA_median_p.scaled`, `linfct_from_model`, `linfct_matrix_contrasts`, `pivot_model_contrasts_2_Wide`, `addContrastResults`, `tidy_FragPipe_combined_protein_deprec`

**Note:** Decorator classes (`LFQDataTransformer`, `LFQDataPlotter`, etc.) are accessed via factory methods — used at runtime but not via `prolfqua::` prefix.

---

## Verification

After any removals:
1. `make document` — regenerate NAMESPACE
2. `make check-fast` — ensure no breakage
3. Grep for removed function names across all `prolfqua_fml/` packages
