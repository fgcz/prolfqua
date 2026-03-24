# TODO: Split tidyMS_R6_Modelling.R (1634 lines)

## Current contents by theme

The file mixes 5 distinct concerns:

### A. Strategy R6 classes + wrappers (lines 1–396)
- `.ehandler` (internal error handler)
- `df.residual.rlm`, `sigma.rlm` (S3 methods for rlm)
- `StrategyLmer` R6 class + `strategy_lmer()` wrapper
- `StrategyLM` R6 class + `strategy_lm()` wrapper
- `StrategyRLM` R6 class + `strategy_rlm()` wrapper
- `AnovaExtractor` R6 class + `get_anova_df()` wrapper

**Consumers:** `build_model()`, `build_model_impute()`, facades, `Model$get_anova()`

### B. Imputation S3 class + fitting internals (lines 399–762)
- `.likelihood_ratio_test` (used by `LR_test` in tidyMS_build_model.R)
- `new_lm_imputed` + S3 methods (`vcov`, `sigma`, `df.residual`)
- `compute_borrowed_variance`
- `impute_refit_singular`
- `isSingular_lm`
- `get_complete_model_fit`
- `model_analyse` (the core per-protein model fitting loop)
- `plot_lmer_peptide_predictions`

**Consumers:** `build_model()`, `build_model_impute()`, `Model$get_anova()`, `Model$get_coefficients()`

### C. Contrast linear function machinery (lines 764–1311)
- `.model_coeff_matrix`, `.get_match_idx`, `.coeff_weights_factor_levels`
- `linfct_from_model`
- `linfct_matrix_contrasts`
- `linfct_all_possible_contrasts`
- `linfct_factors_contrasts`
- `my_contrast` (lm contrast computation)
- `my_contest` (lmer contrast computation via Satterthwaite/KR)
- `contrasts_linfct` (applies contrasts across all proteins)

**Consumers:** `Contrasts` R6 class, `ContrastsROPECA`, strategy `contrast_fun` fields

### D. Moderation + p-value adjustment (lines 1314–1409)
- `moderated_p_limma`
- `moderated_p_limma_long`
- `adjust_p_values`

**Consumers:** `ContrastsModerated`, `Model$get_anova()`, various contrast classes

### E. ROPECA + Fisher exact (lines 1414–1633)
- `get_p_values_pbeta`
- `summary_ROPECA_median_p.scaled`
- `contrasts_fisher_exact`

**Consumers:** `ContrastsROPECA`, `ContrastsMissing`

---

## Suggested split

| New file | Contents | Lines |
|----------|----------|-------|
| `tidyMS_R6_Modelling.R` (keep, shrink) | **A**: Strategy R6 classes + wrappers, `.ehandler`, rlm S3 methods | ~400 |
| `tidyMS_model_fitting.R` (new) | **B**: `model_analyse`, `isSingular_lm`, `get_complete_model_fit`, `impute_refit_singular`, `compute_borrowed_variance`, `new_lm_imputed` + S3, `.likelihood_ratio_test`, `plot_lmer_peptide_predictions` | ~370 |
| `tidyMS_contrasts.R` (new) | **C**: `linfct_*` family, `.model_coeff_matrix`, `my_contrast`, `my_contest`, `contrasts_linfct` | ~550 |
| `tidyMS_moderation.R` (new) | **D + E**: `moderated_p_limma*`, `adjust_p_values`, `get_p_values_pbeta`, `summary_ROPECA_median_p.scaled`, `contrasts_fisher_exact` | ~320 |

This brings each file under 550 lines and groups by clear responsibility:
- Strategies (what model to fit)
- Fitting (how to fit it)
- Contrasts (how to extract comparisons)
- Moderation (how to adjust statistics)

### Alternative: merge B into `tidyMS_build_model.R`

Since `tidyMS_build_model.R` already has `build_model()` and `build_model_impute()` which directly consume `model_analyse`, `impute_refit_singular`, etc., merging B there instead of a new file would co-locate producers and consumers (~600 lines total, still reasonable).
