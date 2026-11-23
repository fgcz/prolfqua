# Fit one model per subject with a modelling strategy

used in project p2901

## Usage

``` r
model_analyse(
  pepIntensity,
  model_strategy,
  subject_id = "protein_Id",
  model_name = "Model",
  label = model_name,
  progress = getOption("prolfqua.progress", NULL)
)
```

## Arguments

- pepIntensity:

  long-format data with one row per observation

- model_strategy:

  modelling strategy, e.g. from
  [`strategy_lm()`](https://fgcz.github.io/prolfqua/reference/strategy.md)

- subject_id:

  columns identifying a subject; one model is fitted per subject

- model_name:

  name stored with the fitted models

- label:

  short label naming this fitting pass in progress output; defaults to
  `model_name`.

- progress:

  progress reporter selector passed to
  [`.make_progress()`](https://fgcz.github.io/prolfqua/reference/dot-make_progress.md);
  `NULL` (default) keeps the legacy terminal
  [`progress::progress_bar`](http://r-lib.github.io/progress/reference/progress_bar.md),
  a `function(i, total, label)` or `"message"` emits batch-visible
  progress. Defaults to `getOption("prolfqua.progress", NULL)`.

## Value

a list with `model_name` and `model_df`, a tibble with one row per
subject holding the nested `data`, the fitted `linear_model` (the error
message for a failed fit), `fit_warnings`, `has_model_fit`,
`isSingular`, `df.residual`, `sigma`, `nr_coef`, `nr_coef_not_NA`,
`converged` and `iterations`

## Note

Deprecated in prolfqua 1.9.5; use
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
