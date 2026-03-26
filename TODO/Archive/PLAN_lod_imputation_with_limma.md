# Plan: LOD Imputation for Limma (`build_model_limma_impute` + `ContrastsLimmaImputeFacade`)

## Context

`build_model_impute()` + `ContrastsLMImputeFacade` handle proteins with missing groups for `lm`: they impute NAs with LOD, refit, and borrow variance from successful fits so that the constant imputation doesn't produce artificially small p-values. Limma has no equivalent. The goal is to give limma the same capability.

The key challenge: limma is matrix-based (`lmFit` fits all proteins at once), so we can't wrap individual models with an S3 class like `lm_imputed`. Instead, we manipulate the `MArrayLM` fit object's per-protein fields directly.

Additionally, we just implemented `ContrastsLimma(eBayes = FALSE)` which returns raw unmoderated statistics — enabling downstream DEqMS moderation on limma results.

## Design Decisions

1. **Only sigma borrowing** (no `borrow_method` parameter). Limma represents per-protein variance as `sigma^2 * stdev.unscaled^2`. The `stdev.unscaled` from the imputed fit captures the design structure; we only need to replace `sigma` with borrowed value. No per-protein vcov matrices exist in MArrayLM.

2. **No new S3 class needed.** MArrayLM is a named list — direct field replacement (`coefficients`, `sigma`, `df.residual`, `stdev.unscaled`) is safe. limma's `contrasts.fit()` and `eBayes()` just read these fields.

3. **New limma-specific helper** `compute_borrowed_variance_limma()` rather than refactoring the existing `compute_borrowed_variance()` which is tightly coupled to the per-protein modelDF tibble.

4. **Degrees of freedom must be corrected for imputed proteins.** The LOD-imputed fit's `fit_lod$df.residual` counts all samples (including imputed ones) as real observations — this is too high and would produce artificially small p-values. We **never** use `fit_lod$df.residual` for imputed proteins. Instead, following `lm_impute`'s `.impute_one_protein()` (line 492-497 of tidyMS_build_model.R):
   - `"observed"`: `df = max(n_observed - p, 1)` where `n_observed = rowSums(!is.na(expr_matrix))` counts only the originally non-missing values per protein, and `p = ncol(design)` is the number of model parameters.
   - `"borrowed"`: `df = median(fit_na$df.residual[good])` from successful proteins.

## Implementation Steps

### Step 1: `compute_borrowed_variance_limma()` in `R/ContrastsLimma.R`

Internal helper (not exported), placed before `build_model_limma()`.

- Input: `fit` (MArrayLM from `lmFit`)
- Identify successful proteins: rows where `!any(is.na(coefficients))`, `df.residual > 1`, `sigma > 0`
- Return `list(sigma = median(fit$sigma[good]), df = median(fit$df.residual[good]))`
- `stop()` if no successful proteins

### Step 2: `build_model_limma_impute()` in `R/ContrastsLimma.R`

Placed after `build_model_limma()`. Signature:

```r
build_model_limma_impute <- function(
  lfqdata, strategy,
  modelName = paste0(strategy$model_name, "Imputed"),
  lod = NULL,
  df_method = c("observed", "borrowed")
)
```

Algorithm:
1. `to_wide(as.matrix = TRUE)` → `expr_matrix`, `annotation`, `rowdata`
2. Build design matrix + resolve weights (reuse same logic as `build_model_limma`)
3. `fit_na <- limma::lmFit(expr_matrix, design, weights = wt)`
4. `failed <- which(rowSums(is.na(fit_na$coefficients)) > 0)`
5. Early return if `length(failed) == 0`
6. Compute LOD if NULL: `MissingHelpers$new(lfqdata$data, lfqdata$config)$get_LOD()`
7. `borrowed <- compute_borrowed_variance_limma(fit_na)`
8. Impute: `expr_imputed <- expr_matrix; expr_imputed[is.na(expr_imputed)] <- lod; expr_imputed <- pmax(expr_imputed, lod)`
9. `fit_lod <- limma::lmFit(expr_imputed, design, weights = wt)`
10. Hybrid fit — replace failed rows in `fit_na`:
    - `fit_na$coefficients[failed, ] <- fit_lod$coefficients[failed, ]`
    - `fit_na$stdev.unscaled[failed, ] <- fit_lod$stdev.unscaled[failed, ]`
    - `fit_na$sigma[failed] <- borrowed$sigma` (borrowed, NOT from `fit_lod`)
    - `fit_na$Amean[failed] <- fit_lod$Amean[failed]` (needed for eBayes trend)
    - **CRITICAL — df correction** (never use `fit_lod$df.residual`, it counts imputed values as real observations):
      - df_method "observed": `fit_na$df.residual[failed] <- pmax(n_observed[failed] - p, 1)` where `n_observed = rowSums(!is.na(expr_matrix))` and `p = ncol(design)`
      - df_method "borrowed": `fit_na$df.residual[failed] <- borrowed$df`
11. Build dummy model (same as `build_model_limma`)
12. Return `ModelLimma$new(fit = fit_na, ...)`

### Step 3: `ContrastsLimmaImputeFacade` in `R/ContrastsFacades.R`

Add after `ContrastsLimmaFacade` (after line 133). Same API surface as all other facades.

```r
ContrastsLimmaImputeFacade <- R6::R6Class(
  "ContrastsLimmaImputeFacade",
  public = list(
    model = NULL,
    contrast = NULL,
    .lfqdata = NULL,
    .contrast_names = NULL,
    initialize = function(lfqdata, modelstr, contrasts,
                          lod = NULL, df_method = c("observed", "borrowed"),
                          weights = lfqdata$config$nr_children, ...) {
      .assert_aggregated_facade_input(lfqdata, "ContrastsLimmaImputeFacade")
      self$.lfqdata <- lfqdata
      self$.contrast_names <- names(contrasts)
      response <- lfqdata$config$get_response()
      full_formula <- paste(response, modelstr)
      strat <- strategy_limma(full_formula, weights = weights, ...)
      self$model <- build_model_limma_impute(lfqdata, strat, lod = lod, df_method = df_method)
      self$contrast <- ContrastsLimma$new(self$model, contrasts)
    },
    get_contrasts = function(...) {
      res <- .add_facade_column(self$contrast$get_contrasts(...), "limma_impute")
      res[!is.na(res$diff), ]
    },
    get_missing = function() .compute_missing(self$.lfqdata, self$.contrast_names, self$get_contrasts()),
    get_Plotter = function(...) self$contrast$get_Plotter(...),
    to_wide = function(...) self$contrast$to_wide(...)
  )
)
```

### Step 4: Register `"limma_impute"` in `R/build_contrast_analysis.R`

- Add `"limma_impute"` to the `method` default vector (line 79)
- Add switch case: `limma_impute = ContrastsLimmaImputeFacade$new(lfqdata, modelstr, contrasts, ...)`

### Step 5: Update `FACADE_REGISTRY` in `R/ContrastsFacades.R`

Add entries (line ~789):
```r
limma_impute = list(class = "ContrastsLimmaImputeFacade", needs = "aggregated"),
lm_impute = list(class = "ContrastsLMImputeFacade", needs = "aggregated")
```

### Step 6: Tests in `tests/testthat/test-ContrastsLimma.R`

1. **`build_model_limma_impute` recovers failed proteins** — simulate with `weight_missing = 0.5`, compare NA coefficient count vs plain `build_model_limma`
2. **`ContrastsLimmaImputeFacade` end-to-end** — `get_contrasts()` returns results, `get_Plotter()` works, `to_wide()` works
3. **Fold changes for non-imputed proteins unchanged** — inner-join on shared proteins, `cor(diff) > 0.99`
4. **`build_contrast_analysis(method = "limma_impute")` dispatches** — check class name
5. **`df_method = "borrowed"` variant** — runs without error
6. **df correction for imputed proteins** — verify that `model$fit$df.residual[failed]` is less than the total sample count minus parameters (i.e., imputed observations are NOT counted as real data points)

### Step 7: `make document`

Regenerate NAMESPACE and man pages.

## Files Modified

| File | Change |
|------|--------|
| `R/ContrastsLimma.R` | Add `compute_borrowed_variance_limma()` + `build_model_limma_impute()` |
| `R/ContrastsFacades.R` | Add `ContrastsLimmaImputeFacade`; update `FACADE_REGISTRY` |
| `R/build_contrast_analysis.R` | Add `"limma_impute"` to switch + method arg |
| `tests/testthat/test-ContrastsLimma.R` | Add 5 test cases |

## Verification

```bash
make document                           # regenerate NAMESPACE
make install                            # install with new exports
make test                               # full test suite
Rscript -e "testthat::test_file('tests/testthat/test-ContrastsLimma.R')"  # focused
```

Smoke test in R:
```r
library(prolfqua)
istar <- sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5)
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

# Plain limma (drops proteins with missing groups)
fa1 <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma")
nrow(fa1$get_contrasts())

# Limma with LOD imputation (recovers them)
fa2 <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma_impute")
nrow(fa2$get_contrasts())  # should be >= fa1

# Fold changes for shared proteins should correlate highly
merged <- dplyr::inner_join(
  dplyr::select(fa1$get_contrasts(), protein_Id, diff1 = diff),
  dplyr::select(fa2$get_contrasts(), protein_Id, diff2 = diff),
  by = "protein_Id"
)
cor(merged$diff1, merged$diff2, use = "complete.obs")  # > 0.99
```
