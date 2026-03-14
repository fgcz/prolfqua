---
name: adding-models-to-prolfqua
description: Add or review modelling backends in prolfqua. Use when implementing a new model adapter, contrast adapter, or facade, and when checking compatibility with ModelInterface and ContrastsInterface.
---
# Adding Models To Prolfqua
Use this skill when extending the modelling API in `prolfqua`.

The main goal is not choosing an internal fitting style. The main goal is implementing the public contracts correctly:
- `ModelInterface`
- `ContrastsInterface`
- facade classes in `R/ContrastsFacades.R`

## Start With The Public API
A backend is integrated successfully only when all three layers are coherent.

1. A model builder creates the fitted backend object.
2. A `Model*` class adapts that object to `ModelInterface`.
3. A `Contrasts*` class adapts hypothesis testing results to `ContrastsInterface`.
4. A facade provides the user-facing workflow.

## Step 1: Prepare LFQData Correctly
Before adding code, verify:
- the response column: `lfqdata$config$get_response()`
- the analysis unit: `lfqdata$subject_Id()`
- required factors exist in `lfqdata$factors()`
- the data is already transformed if the backend expects log scale

Only then decide whether the backend consumes nested subject-wise data or `lfqdata$to_wide(as.matrix = TRUE)`. That choice is secondary to the interface contract.

## Step 2: Implement The Model Adapter First
Every new backend should have a `Model*` adapter inheriting `ModelInterface`.

Required methods:
- `get_coefficients()`
- `get_anova()`
- `coef_histogram()`
- `coef_volcano()`
- `coef_pairs()`
- `anova_histogram()`

Expected `get_coefficients()` shape:
- subject ID columns
- `factor`
- effect column, usually `Estimate`
- p-value column for plotting, usually `Pr...t..`

Expected `get_anova()` shape:
- subject ID columns
- `factor`
- `p.value`
- `FDR`

Rules:
- Keep column names compatible with the existing `Model` API.
- If the backend only supports an omnibus test, state that explicitly in docs and examples.
- Do not leak backend-specific result shapes into downstream code.

## Step 3: Implement The Contrast Adapter
Every new backend should have a `Contrasts*` adapter inheriting `ContrastsInterface`.

Required methods:
- `get_contrast_sides()`
- `get_contrasts()`
- `get_Plotter()`
- `to_wide()`

`get_contrasts()` should return this standard schema:
- `modelName`
- `contrast`
- `avgAbd`
- `diff`
- `FDR`
- `statistic`
- `std.error`
- `df`
- `p.value`
- `conf.low`
- `conf.high`
- `sigma`

Rules:
- Translate backend-specific names inside the adapter.
- Keep downstream code unaware of backend-specific output conventions.
- Reuse `pivot_model_contrasts_2_Wide()` for `to_wide()` unless a backend truly requires something else.
- Validate the final output against `ContrastsInterface$column_description()`.

## Step 4: Only Then Decide How To Fit The Model
After the interface contract is clear, choose the fitting implementation.

### Reuse the classic path
If the backend can reuse `build_model()`, add a `strategy_*()` like `strategy_lm()` or `strategy_lmer()`.

The strategy list must contain:
- `model_fun`
- `isSingular`
- `contrast_fun`
- `model_name`
- `report_columns`
- `anova_df`
- `is_mixed`
- `df_residual`
- `sigma`

Template:
```r
model_fun <- function(x, pb, get_formula = FALSE) {
  if (get_formula) return(formula)
  if (!missing(pb)) pb$tick()
  tryCatch(fit_backend(formula, data = x), error = .ehandler)
}
```

Rules:
- failures must return a character string via `.ehandler`, not `NULL`
- `model_analyse()` uses `!is.character(x)` as success
- `contrast_fun` must return tidy contrast statistics

### Add a dedicated builder
If the backend does not fit `model_analyse()`, create:
- `build_model_<backend>()`
- `Model<Backend>`
- `Contrasts<Backend>`

Use this for wide-matrix, Bayesian, or backend-specific pipelines.

## Step 5: Make Contrast Construction Reliable
This is the main place integrations break.

The default `Contrasts` path assumes:
1. `linfct_from_model()` can recover coefficient structure
2. `linfct_matrix_contrasts()` can parse user expressions like `"group_A - group_B"`
3. the backend can evaluate those linear functions correctly

Preferred rule:
- derive contrast structure from the backend's design matrix or coefficient names
- do not infer the global contrast scaffold from one arbitrary incomplete row

If those assumptions fail, do not force the backend into `Contrasts`. Write a dedicated `Contrasts*` adapter.

## Step 6: Add The Facade
If the backend should be part of the user-facing API, add a facade in `R/ContrastsFacades.R`.

Facade responsibilities:
1. validate the shape of `LFQData`
2. prepend the response column to the formula
3. fit the model
4. build the contrast object
5. return standardized output with a `facade` column

Follow existing names like:
- `ContrastsLMFacade`
- `ContrastsLmerFacade`
- `ContrastsLimmaFacade`

The facade is the final integration point. If the facade feels awkward, the underlying adapter design is usually still wrong.

## Step 7: Test Interface Compliance
At minimum, add tests for:

1. construction
- the builder or strategy returns the expected type

2. `ModelInterface`
- `get_coefficients()` returns documented columns
- `get_anova()` returns documented columns

3. `ContrastsInterface`
- `get_contrasts()` returns the standard schema
- `to_wide()` works
- `get_Plotter()` works

4. edge cases
- missingness
- rank deficiency
- multi-factor designs
- multiple contrasts

5. invariants
- fold changes match a trusted implementation when possible
- output schema is stable across backends

## Step 8: Validate The Common Failure Modes
Explicitly test:
- sample order alignment for wide backends
- correct grouping unit for subject-wise backends
- dropped coefficients and absent levels
- interaction terms
- final output schema and column names

## Step 9: Document And Export
When the backend is ready:
1. add roxygen docs in `R/*.R`
2. add examples showing the full pipeline
3. run `make document`
4. never edit `NAMESPACE` directly

## Review Checklist
- [ ] `ModelInterface` is fully implemented
- [ ] `ContrastsInterface` is fully implemented
- [ ] facade API is coherent and easy to call
- [ ] backend-specific columns are normalized inside adapters
- [ ] contrast construction is stable
- [ ] tests cover missingness, rank deficiency, and multi-factor designs
- [ ] `make document` was run after roxygen edits

## Files To Study First
- `R/ModelInterface.R`
- `R/ContrastsInterface.R`
- `R/Model.R`
- `R/Contrasts.R`
- `R/ContrastsLimma.R`
- `R/ContrastsFacades.R`
- `R/LFQData.R`
- `R/tidyMS_R6Model.R`
- `R/tidyMS_R6_Modelling.R`
- `tests/testthat/test-Model.R`
- `tests/testthat/test-Contrasts.R`
- `tests/testthat/test-ContrastsLimma.R`
