# LM-based Imputation with Covariance Borrowing — Implementation Design

## Problem

When a protein has too many missing values, `lm()` either:
1. **Fails entirely** — returns a character error string (`exists_lmer = FALSE`)
2. **Produces NA coefficients** — a group has some but insufficient data (`isSingular = TRUE`)
3. **Drops a group entirely** — a group has zero observations, so the corresponding coefficient simply doesn't exist in the model (`nrcoef < max(nrcoef)`)

In all three cases, contrasts involving the affected groups cannot be estimated. The existing fallback (`ContrastsMissing`) uses a completely separate group-mean approach that bypasses the lm infrastructure. We wanted an alternative that stays within the lm framework.

## Core Idea

```
For each protein:
  1. Try normal lm fit
  2. If success and complete → keep fit, contribute to "donor pool"
  3. If failure/singular/incomplete coefficients:
     a. Complete data with all sample rows (add missing group rows)
     b. Replace NAs with LOD, clamp all values to max(value, LOD)
     c. Fit lm on imputed data → get correct point estimates (coefficients)
     d. Replace vcov/sigma/df with borrowed values from donor pool
  4. Extract contrasts using standard Contrasts machinery (unchanged)
```

The imputation gives correct fold-change *direction* — the LOD-clamping ensures a group with all imputed values can never appear higher than a group with observed values. The borrowed covariance gives honest *uncertainty* — since we imputed with a constant, the model's own residual variance is artificially small, so we substitute a typical variance from proteins that fit successfully.

## Architecture

### Key Insight: S3 Dispatch Trick

The entire contrast computation pipeline (`my_contrast_V2` → `my_contrast`) calls three S3 generics on the model object:
- `vcov(m)` — to get the variance-covariance matrix
- `sigma(m)` — to get the residual standard error
- `df.residual(m)` — to get the residual degrees of freedom

By wrapping the refitted lm in an `lm_imputed` S3 class that overrides these three methods, the existing contrast code works *completely unchanged*. `coefficients(m)` still dispatches to the underlying `lm` class (inherited), returning the coefficients from the imputed fit.

### Separation of Concerns

- `build_model()` — unchanged, the original function
- `build_model_impute()` — new, dedicated function for imputed models
- `Contrasts`, `ContrastsModerated`, `my_contrast_V2` — all unchanged, completely unaware of imputation

## Components

### 1. `lm_imputed` S3 Class (`R/tidyMS_R6_Modelling.R`)

**Constructor:** `new_lm_imputed(model, borrowed_vcov, borrowed_sigma, borrowed_df, n_observed)`

Stores borrowed values as attributes on the lm object and prepends `"lm_imputed"` to the class vector. The three S3 method overrides return the borrowed values:

```r
vcov.lm_imputed     → attr(object, "borrowed_vcov")
sigma.lm_imputed    → attr(object, "borrowed_sigma")
df.residual.lm_imputed → attr(object, "borrowed_df")
```

All three methods are registered in NAMESPACE via `#' @export` roxygen tags (roxygen detects the `generic.class` pattern).

### 2. `compute_borrowed_variance(modelDF, method)` (`R/tidyMS_R6_Modelling.R`)

Computes the donor pool statistics from successful, non-singular, complete models in `modelDF`.

**Two methods (user-selectable):**

- `method = "sigma"` (Option A): Returns `median(sigma)` and `median(df.residual)` across donors. The actual vcov for each imputed protein is computed later as `borrowed_sigma² × (X'X)⁻¹` using the imputed model's own design matrix. This preserves the per-protein design structure while borrowing the variance scale.

- `method = "vcov"` (Option B): Stacks all donor vcov matrices into a 3D array and takes element-wise median. Falls back to sigma method if donors have different coefficient dimensions (with a warning).

**Donor selection:** `get_complete_model_fit(modelDF)` filtered to `isSingular == FALSE`. This selects proteins with `exists_lmer == TRUE`, full coefficient rank (`nrcoeff_not_NA == max`), and `df.residual > 1`.

### 3. `impute_refit_singular(...)` (`R/tidyMS_R6_Modelling.R`)

The core loop. Identifies which proteins need imputation, then processes each one.

**Detection — three conditions (OR):**
```r
needs_impute <- (!modelDF$exists_lmer) |                          # lm() errored
  (!is.na(modelDF$isSingular) & modelDF$isSingular) |             # NA coefficients or df < 2
  (!is.na(modelDF$nrcoef) & modelDF$nrcoef < max_coef)           # missing group(s)
```

The third condition was the non-obvious one. When an entire group is absent from a protein's data, `lm()` fits successfully with fewer coefficients (e.g., 2 instead of 3 for a 3-group design). `isSingular` is FALSE because there are no NA coefficients — the coefficient simply doesn't exist. But contrasts referencing the missing group fail silently (produce NA results in `my_contrast_V2`).

**Per-protein steps:**

1. **Complete data with sample template** — `left_join(sample_template, dat)` adds rows for all samples, filling the response with NA for samples not present. This is critical for the "missing group" case where the nested data has zero rows for one or more groups.

2. **Impute** — `ifelse(is.na(response), lod, response)` fills NAs with LOD.

3. **Clamp** — `pmax(response, lod)` ensures no observed value is below LOD. This prevents the situation where an observed group has a lower estimate than an imputed group (which would produce a wrong-direction fold change).

4. **Refit** — `model_strategy$model_fun(dat)` fits a new lm on the complete imputed data. Skip if it still fails.

5. **Compute borrowed vcov** — For sigma method: `borrowed_sigma² × summary(new_model)$cov.unscaled`. For vcov method: use the pre-computed element-wise median.

6. **Compute df** — Two options:
   - `df_method = "observed"`: `max(n_observed - p, 1)` — honest but can be very small
   - `df_method = "borrowed"`: median df from donors — more stable, similar to limma's prior df

7. **Wrap** — `new_lm_imputed(new_model, ...)` creates the S3 wrapper.

8. **Update modelDF** — Replace `linear_model`, `data`, `exists_lmer`, `isSingular`, `sigma`, `df.residual`, `nrcoef`, `nrcoeff_not_NA` for the protein.

### 4. `build_model_impute(lfqdata, model_strategy, ...)` (`R/tidyMS_R6Model.R`)

Dedicated entry point. Requires an `LFQData` object (not a raw data.frame) because it needs:
- `lfqdata$config$get_response()` for the response column name
- `lfqdata$data` to build the sample template (all sample × group combinations)
- `MissingHelpers$new(...)$get_LOD()` for automatic LOD computation

**Sample template construction:**
```r
non_subject_cols <- setdiff(colnames(lfqdata$data), c(subject_Id, response))
sample_template <- lfqdata$data |> select(all_of(non_subject_cols)) |> distinct()
```
This gives all unique combinations of sample-level columns (sampleName, group_, etc.) excluding the protein identifier and response. When left-joined into a protein's nested data, it adds rows for every sample, with NA response for missing ones.

**LOD auto-computation:** Delegates to `MissingHelpers$get_LOD()`, which takes the median of group means where `nrMeasured == 1` (groups with exactly one observation — likely just above the detection limit).

### 5. `ContrastsLMImputeFacade` (`R/ContrastsFacades.R`)

Follows the exact same pattern as `ContrastsLMFacade`:

```r
initialize:
  strategy_lm(formula)  →  build_model_impute(lfqdata, strat)  →
  ContrastsModerated$new(Contrasts$new(model, contrasts))
```

Exposes: `get_contrasts()`, `get_missing()`, `get_Plotter()`, `to_wide()`.

Registered in `build_contrast_analysis()` as `method = "lm_impute"`.

### 6. Integration into `build_contrast_analysis()` (`R/build_contrast_analysis.R`)

Added `"lm_impute"` to the method enum and switch:
```r
lm_impute = ContrastsLMImputeFacade$new(lfqdata, modelstr, contrasts, ...)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `lod` | NULL (auto) | Limit of detection. If NULL, computed from data via MissingHelpers. |
| `borrow_method` | `"sigma"` | `"sigma"`: borrow scalar σ², use per-protein (X'X)⁻¹. `"vcov"`: borrow element-wise median of full vcov matrices. |
| `df_method` | `"observed"` | `"observed"`: max(n_observed - p, 1). `"borrowed"`: median df from successful fits. |

## Files Modified

| File | What |
|------|------|
| `R/tidyMS_R6_Modelling.R` | `new_lm_imputed`, 3 S3 methods, `compute_borrowed_variance`, `impute_refit_singular` |
| `R/tidyMS_R6Model.R` | `build_model_impute` (new function; `build_model` unchanged) |
| `R/ContrastsFacades.R` | `ContrastsLMImputeFacade` |
| `R/build_contrast_analysis.R` | Added `"lm_impute"` method |
| `vignettes/ContrastFacades.Rmd` | Added `lm_impute` to the vignette |
| `tests/testthat/test-ImputeModel.R` | 6 tests, 20 assertions |

**Unchanged:** `Contrasts.R`, `ContrastsModerated.R`, `Model.R`, `ContrastsInterface.R`, `my_contrast_V2`, `my_contrast`, `contrasts_linfct` — the entire contrast computation pipeline is unaware of imputation.

## Data Flow Diagram

```
LFQData (protein-level, has missing values)
    │
    ▼
build_model_impute(lfqdata, strategy_lm(...))
    │
    ├─ model_analyse()                    # First pass: fit all proteins
    │   └─ modelDF with:
    │       ├─ 77 complete models (exists_lmer=T, isSingular=F, nrcoef=3)
    │       └─  3 incomplete models (nrcoef=2, missing group)
    │
    ├─ MissingHelpers$get_LOD()           # Compute LOD from data
    │
    ├─ Build sample_template              # All 12 sample rows
    │
    └─ impute_refit_singular()            # Second pass: fix the 3
        │
        ├─ compute_borrowed_variance()    # median(sigma) from 77 donors
        │
        └─ For each of the 3 proteins:
            ├─ left_join(sample_template)  # Add 4 missing sample rows
            ├─ Impute NAs with LOD         # Fill the 4 new rows
            ├─ Clamp to max(value, LOD)    # Floor all values
            ├─ Refit lm()                  # Now has 3 coefficients
            ├─ Compute borrowed vcov       # σ²_borrowed × (X'X)⁻¹
            └─ Wrap in lm_imputed          # S3 override for vcov/sigma/df
    │
    ▼
Model (80 proteins, all with full coefficients)
    │
    ▼
Contrasts$new(model, contrasts)           # Standard, unchanged
    │
    ▼
ContrastsModerated$new(contrasts)         # Standard, unchanged
    │
    ▼
160/160 contrast results (80 proteins × 2 contrasts, 0 missing)
```

## Lessons Learned

1. **Missing data has three faces** — not just NA coefficients (`isSingular`) or failed fits (`exists_lmer`), but also entirely absent groups where `lm()` succeeds with fewer coefficients. The `nrcoef < max(nrcoef)` check was the non-obvious third condition.

2. **Nested data is already filtered** — `model_analyse()` nests data by protein *after* NAs are removed. A protein missing an entire group has zero rows for that group in its nested data. Simply replacing NAs doesn't help — you need to add the missing rows first via a sample template join.

3. **S3 dispatch is the cleanest injection point** — rather than modifying contrast functions, passing extra parameters through the strategy, or creating wrapper classes, overriding `vcov()`/`sigma()`/`df.residual()` keeps all changes in one place and the rest of the pipeline untouched.

4. **Keep `build_model` clean** — rather than adding `impute = FALSE` with conditional logic, a separate `build_model_impute` function is simpler and more explicit.
