# Summarize modelling and error reporting

Summarize modelling and error reporting

## Usage

``` r
model_summary(mod)
```

## Arguments

- mod:

  model table see
  [`build_model`](https://fgcz.github.io/prolfqua/reference/build_model.md)

## Value

a list with the tables `exists` (counts of subjects with and without a
model fit) and `isSingular` (counts of singular and non-singular fits)

## Note

Deprecated in prolfqua 1.9.5; use
[`build_contrast_analysis()`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md).
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
