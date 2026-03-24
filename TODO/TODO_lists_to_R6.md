# TODO: Convert list-based pseudo-classes to R6

## Background

Several functions in the codebase construct and return named lists that act as
pseudo-objects — carrying functions (closures), config parameters, or structured
return values. These are candidates for conversion to proper R6 classes to
improve discoverability, type safety, and consistency with the rest of the
codebase.

The candidates are organized by priority: **Tier 1** are true "list-based
classes" (contain functions or act as behavioural objects), **Tier 2** are
structured return-value bundles that would benefit from a lightweight R6 wrapper.

---

## Tier 1 — List-based classes (contain functions / closures)

### 1. `strategy_limma()` → `StrategyLimma` R6 class

- **File:** `R/ContrastsLimma.R:21`
- **Current shape:** Returns `list(formula, model_name, trend, robust, weights)`
- **Why convert:** Every other strategy is already an R6 class (`StrategyLM`,
  `StrategyRLM`, `StrategyLmer`, `StrategyLogistf`). `strategy_limma()` is the
  only one that returns a plain list, creating an inconsistency in the strategy
  interface. `build_model_limma()` accesses members via `$` on this list.
- **Suggested R6 fields:** `formula`, `model_name`, `trend`, `robust`, `weights`
- **Consumers:** `build_model_limma()`, `ContrastsLimmaFacade`, `ContrastsDEqMSFacade`
- **Complexity:** Low — pure data, no methods needed initially. Can add
  `model_fun`, `isSingular`, etc. later to unify with other Strategy classes.

### 2. `get_anova_df()` → `AnovaExtractor` R6 class (or similar)

- **File:** `R/tidyMS_R6_Modelling.R:344`
- **Current shape:** Returns `list(fun = <closure>, col_pval = "Pr..F..", col_fdr = "FDR.Pr..F..")`
- **Why convert:** This is a classic "function + metadata" bundle. The closure
  `fun` takes a model and returns an ANOVA data frame; `col_pval` and `col_fdr`
  describe the output columns. Wrapping in R6 makes the contract explicit and
  allows adding methods (e.g., `$extract(model)` instead of `$fun(model)`).
- **Suggested R6 fields:** `test` (stored param)
- **Suggested R6 methods:** `$extract(model)` (currently `$fun`), `$col_pval()`,
  `$col_fdr()`
- **Consumers:** Used in model summary pipelines, `Model` class
- **Complexity:** Low

### 3. `hierarchy_counts_sample()` → `HierarchyCountsSample` R6 class

- **File:** `R/AnalysisConfiguration.R:748`
- **Current shape:** Returns a single closure `function(value = c("wide", "long", "plot"))`
  that captures `summary` (a data frame) and `configuration`.
- **Why convert:** The closure dispatches on a string argument to return
  different views (wide df, long df, ggplot). An R6 class with `$wide()`,
  `$long()`, `$plot()` methods would be clearer and more discoverable.
- **Suggested R6 fields:** `summary` (data frame), `configuration`
- **Suggested R6 methods:** `$wide()`, `$long()`, `$plot()`
- **Consumers:** `LFQDataSummariser`
- **Complexity:** Low

---

## Tier 2 — Structured return bundles (data + metadata)

These functions return `list(...)` as a way to return multiple values. They
don't contain closures, but wrapping them in lightweight R6 (or S3) classes
would add type clarity and prevent field-name typos.

### 4. `model_analyse()` return value → `ModelResult`

- **File:** `R/tidyMS_R6_Modelling.R:641`
- **Returns:** `list(modelDF = <tibble>, modelName = <string>)`
- **Consumers:** `build_model()`, `build_model_impute()`, `Model$new()`,
  `model_summary()`, downstream contrasts code
- **Note:** This is the most widely consumed return-value bundle. Already has an
  implicit contract (`$modelDF`, `$modelName`). The `Model` R6 class wraps this,
  so converting the raw return might be redundant — but making `Model$new()`
  accept either form would be a nice transitional step.

### 5. `compute_borrowed_variance()` return value → `BorrowedVariance`

- **File:** `R/tidyMS_R6_Modelling.R:458`
- **Returns:** `list(vcov?, sigma, df, method)` (shape depends on method)
- **Consumers:** `build_model_impute()`, `new_lm_imputed()`
- **Note:** The conditional shape (vcov present only when `method == "vcov"`)
  is a good reason to formalize — an R6 class could have `$has_vcov()`.

### 6. `summarize_stats_quantiles()` return value

- **File:** `R/tidyMS_stats.R:290`
- **Returns:** `list(long = <tibble>, wide = <tibble>)`
- **Consumers:** Sample size estimation, QC reporting
- **Complexity:** Very low, borderline — might not be worth R6 overhead.

### 7. `lfq_power_t_test_quantiles()` return value

- **File:** `R/tidyMS_stats.R:395`
- **Returns:** `list(long = <tibble>, summary = <tibble>)`
- **Consumers:** Power/sample-size reporting
- **Complexity:** Same as above — borderline.

### 8. `.add_nr_children()` return value

- **File:** `R/tidyMS_aggregation.R:544`
- **Returns:** `list(data = <tibble>, config = <AnalysisConfiguration>)`
- **Consumers:** `LFQDataAggregator`
- **Note:** Internal function. Conversion optional.

### 9. `.model_coeff_matrix()` return value

- **File:** `R/tidyMS_R6_Modelling.R:726`
- **Returns:** `list(mm = <matrix>, coeffs = <named numeric>)`
- **Consumers:** `linfct_from_model()`, internal contrast generation
- **Note:** Internal function. Conversion optional.

### 10. Missingness functions return values

- **Files:** `R/tidyMS_missigness.R:187, 224, 256, 286`
- **Returns:** Various `list(data = ..., figure = ...)` or `list(data = ..., nsets = ...)`
- **Consumers:** `MissingHelpers`, `LFQDataSummariser`
- **Note:** Several related functions with similar shapes — could share a
  common `MissingnessResult` class.

---

## Recommended order of conversion

1. **`strategy_limma()` → `StrategyLimma`** — highest impact, aligns with
   existing Strategy R6 pattern, low effort
2. **`get_anova_df()`** — small, self-contained, good practice target
3. **`hierarchy_counts_sample()`** — small, improves API clarity
4. Tier 2 items as needed during related refactoring

## Notes

- Keep backward-compatible wrapper functions (e.g., `strategy_limma()` returns
  `StrategyLimma$new(...)`) so existing call sites don't break.
- For `StrategyLimma`, consider whether it should eventually implement the same
  method interface as `StrategyLM` (`model_fun`, `isSingular`, `contrast_fun`,
  `df_residual`, `sigma`) to fully unify the strategy pattern. This is a bigger
  design decision since `build_model_limma()` works fundamentally differently
  from `model_analyse()`.
