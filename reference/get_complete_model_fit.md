# Retrieve complete model fits

Retrieve complete model fits

## Usage

``` r
get_complete_model_fit(complete_models)
```

## Arguments

- complete_models:

  model table, the `model_df` of
  [`model_analyse()`](https://fgcz.github.io/prolfqua/reference/model_analyse.md)

## Value

the rows of `complete_models` that have a model fit, the maximal number
of estimated coefficients and more than one residual degree of freedom

## Note

Deprecated in prolfqua 1.9.5; use
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
