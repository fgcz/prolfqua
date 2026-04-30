# TODO: Long Functions

Generated from the current `prolfqua/R` sources with BiocCheck 1.46.3 function-length logic
(`BiocCheck:::getFunctionLengths`). The current sources contain 37 functions longer than 50 lines.

Note: `TODO_BiocCheck.md` currently says the longest function is `fitFDistRobustly_LG()` at 267 lines. In the current
source, the longest function is `fitFDist_LG()` at 136 lines; `fitFDistRobustly_LG()` has already been shortened below
the 50-line threshold.

## Should Refactor

These functions are either central workflows, contain duplicated logic, or have enough branching that helper extraction
would reduce maintenance risk.

| Function | Lines | File | Start | End | Plan to shorten |
| --- | ---: | --- | ---: | ---: | --- |
| `setup_analysis` | 128 | `R/tidyMS_data_setup.R` | 51 | 178 | Extract hierarchy unification, factor unification, sample-name creation, implicit-column defaults, duplicate-observation check, and complete-case handling into private helpers. Add focused tests before changing because this is central setup code. |
| `build_model_limma_voom_impute` | 127 | `R/ContrastsLimma.R` | 513 | 639 | Promote the nested `.voom_fit` helper to a shared private helper used by both voom builders. Reuse a common limma-imputation helper for failed-row detection, LOD filling, borrowed variance, df correction, dummy-model creation, and `ModelLimma$new()` assembly. |
| `get_contrasts` | 103 | `R/ContrastsLimma.R` | 951 | 1053 | Split the R6 method into helpers for contrast-matrix construction, `limma::contrasts.fit`/`eBayes`, per-contrast result extraction, raw-statistic extraction, average-abundance calculation, p-value adjustment, and cache assignment. |
| `build_model_limma_impute` | 92 | `R/ContrastsLimma.R` | 289 | 380 | Share the failed-fit hybridization workflow with `build_model_limma_voom_impute`: compute LOD, impute matrix/weights, refit, replace failed rows, apply df policy, and construct the dummy model in private helpers. |
| `get_contrasts` | 85 | `R/ContrastsModeratedDEqMS.R` | 340 | 424 | Extract count joining/sanitization, DEqMS moderation, compact-vs-full column shaping, p-value adjustment, and model-name suffixing into helpers shared with `ContrastsModerated`. |
| `.elist_to_lfqdata` | 82 | `R/LFQDataAggregator.R` | 312 | 393 | Split row-id parsing, matrix-to-long expansion, annotation join-column selection, reduced-config creation, and `nr_children`/SE metadata assignment into private helpers. |
| `linfct_matrix_contrasts_vectorized` | 80 | `R/tidyMS_contrasts.R` | 218 | 297 | Extract contrast-name normalization, expression parsing, fast mutate path, fallback mutate path, and failure-warning formatting. Consider sharing fallback warning code with `linfct_matrix_contrasts`. |
| `.rlm_estimate` | 78 | `R/tidyMS_aggregation.R` | 331 | 408 | Extract degenerate single-feature handling, single-sample handling, singular model-matrix reduction, robust fit summarization, invalid-weight cleanup, and final sample join. |
| `moderated_p_deqms_long` | 75 | `R/ContrastsModeratedDEqMS.R` | 150 | 224 | Extract per-group moderation with warning capture, group-label construction, warning-summary formatting, and final bind. This will make the main function a group-map wrapper. |
| `compute_contrast_vectorized` | 72 | `R/tidyMS_contrasts.R` | 487 | 558 | Split empty-result construction, coefficient/vcov alignment, NA-coefficient invalid-row detection, inferential-statistic calculation, and output table assembly. |
| `make_contrasts` | 71 | `R/tidyMS_contrasts.R` | 221 | 291 | This is a nested helper inside `linfct_matrix_contrasts_vectorized`; move it to a private top-level helper and split fast-path evaluation from fallback evaluation and warning reporting. |
| `build_model_limma_voom` | 67 | `R/ContrastsLimma.R` | 415 | 481 | Reuse the same private voom-weight helper planned for `build_model_limma_voom_impute`. Extract diagnostic plotting and `ModelLimma$new()` construction. |
| `model_analyse` | 64 | `R/tidyMS_build_model.R` | 615 | 678 | Extract model fitting, fit-status annotation, coefficient counters, and model-summary column calculation. Replace local helper functions with private top-level helpers if reused. |
| `build_model_impute` | 60 | `R/tidyMS_build_model.R` | 180 | 239 | Extract sample-template construction, `nr_children` resolution, LOD resolution, and final `Model$new()` assembly. Share helper names with `.impute_one_protein` where practical. |
| `get_contrasts` | 59 | `R/ContrastsModerated.R` | 94 | 152 | Extract moderation call, compact/full output shaping, p-value adjustment, and model-name suffixing. Share the column-shaping helper with `ContrastsModeratedDEqMS$get_contrasts()`. |
| `impute_refit_singular` | 58 | `R/tidyMS_build_model.R` | 421 | 478 | Extract singular-row selection, borrowed-variance computation, per-protein imputation mapping, successful-result filtering, and model-table replacement. |
| `contrasts_linfct` | 58 | `R/tidyMS_contrasts.R` | 654 | 711 | Extract linfct normalization, per-model contrast computation with progress, failed-contrast warning/filtering, unnesting, and model-info join. |
| `linfct_matrix_contrasts` | 57 | `R/tidyMS_contrasts.R` | 152 | 208 | Move nested `make_contrasts` to a shared helper, reuse contrast-name/failure-warning helpers from the vectorized path, and leave the wrapper to transpose in/out. |
| `aggregate` | 52 | `R/LFQDataAggregator.R` | 463 | 514 | Extract DPC estimation, impute-only `dpcQuantByRow`, aggregation-mode protein-id construction, aggregation-mode `dpcQuant`, rowname restoration, and LFQData conversion. |

## Nice To Refactor

These are understandable today, but helper extraction would improve readability, test focus, or reuse.

| Function | Lines | File | Start | End | Plan to shorten |
| --- | ---: | --- | ---: | ---: | --- |
| `sim_lfq_data` | 96 | `R/simulate_LFQ_data.R` | 28 | 123 | Separate protein table simulation, fold-change assignment, optional peptide expansion, long-format conversion, replicate sampling, and final sample labeling. Keep the exported function as a readable pipeline. |
| `moderated_p_deqms` | 76 | `R/ContrastsModeratedDEqMS.R` | 47 | 122 | Split finite-row validation, LOESS fit/prediction, prior-df estimation, posterior-variance calculation, moderated-statistic calculation, and result-column attachment. |
| `plot_hierarchies_boxplot` | 72 | `R/tidyMS_plotting.R` | 141 | 212 | Extract plot-data preparation, base plot setup, color/scale handling, beeswarm layer construction, mean-summary layers, and optional faceting. |
| `sim_build_models_logistf` | 67 | `R/logistf.R` | 229 | 295 | Extract simulation-data selection for protein/peptide cases, binary-response encoding, model-formula selection, and final `build_model_logistf()` call. |
| `plot_hierarchies_line_default` | 63 | `R/tidyMS_aggregation.R` | 3 | 65 | Extract facet-formula construction, isotope vs non-isotope aesthetics, point/line layer creation, and final axis/title/log-scale styling. |
| `.impute_one_protein` | 63 | `R/tidyMS_build_model.R` | 494 | 556 | Split sample-template completion, `nr_children` filling, LOD imputation/clamping, refit/error handling, df policy, borrowed covariance, and return-row assembly. |
| `.volcano` | 63 | `R/utilities.R` | 359 | 421 | Extract vertical/horizontal shape builders, base `plot_ly()` creation, threshold-shape list construction, and annotation layout. |
| `plot_pca` | 57 | `R/tidyMS_plotting.R` | 538 | 594 | Extract missing-feature filtering, PCA score/variance calculation, PC bounds check, annotation join, aesthetic/layer selection, and shape-scale addition. |
| `contrasts_linfct_firth` | 56 | `R/logistf.R` | 20 | 75 | Align structure with `contrasts_linfct`: extract per-model contrast computation, failed-result reporting/filtering, unnesting, and model-info join. |
| `contrasts_fisher_exact` | 56 | `R/tidyMS_moderation.R` | 272 | 327 | Move relative-risk, odds-ratio, and single-row Fisher test helpers to private top-level functions. Replace the explicit loop with `purrr::pmap()` or a small helper returning the bound result. |
| `sim_lfq_data_2factor_config` | 55 | `R/simulate_LFQ_data.R` | 281 | 335 | Extract the nested `simulate` body into helpers for base simulation, factor assignment, missingness application, config construction, and seeded execution. |
| `get_contrast_estimates` | 54 | `R/tidyMS_missingness_imputation.R` | 134 | 187 | Extract imputed-abundance contrast creation, missing-indicator contrast creation, measured-count contrast creation, and direction-based estimate correction. |
| `.multigroup_volcano` | 54 | `R/utilities.R` | 172 | 225 | Extract base ggplot construction, optional color/text aesthetics, threshold layers, and faceting/labelling. Also handle the palette setup separately. |
| `multigroup_volcano` | 51 | `R/utilities.R` | 120 | 170 | Extract label preparation, significant-label subset selection, and text-repel layer construction. Keep the exported function as a thin wrapper around `.multigroup_volcano()`. |

## Good As Is

These are low-priority despite the line count: either copied numerical routines where refactoring adds risk, or small
linear functions only just above the threshold.

| Function | Lines | File | Start | End | Plan to shorten |
| --- | ---: | --- | ---: | ---: | --- |
| `fitFDist_LG` | 136 | `R/squeezeVarRob.R` | 64 | 199 | Leave unless BiocCheck compliance becomes mandatory. This is a limma-derived numerical routine; splitting it may obscure the reference algorithm and increase regression risk. |
| `trigammaInverse` | 60 | `R/squeezeVarRob.R` | 3 | 62 | Leave unless changing the numerical code for another reason. The guard clauses and Newton iteration are compact enough and easier to audit in one place. |
| `build_model_limpa` | 51 | `R/ContrastsLimpa.R` | 104 | 154 | Leave for now. It is linear setup/validation/fitting code, only one line over the threshold, and extraction would add more indirection than clarity. |
| `transform_work_intensity` | 51 | `R/LFQDataTransformer.R` | 261 | 311 | Leave for now. It is only one line over the threshold and already has a simple flow: name resolution, diagnostics, mutate, return. |
