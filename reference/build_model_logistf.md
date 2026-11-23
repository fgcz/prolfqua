# Fit Firth logistic regression to binary detection data

Fits Firth's logistic regression per subject to the binary detection
response of aggregated data, one row per subject and sample.

## Usage

``` r
build_model_logistf(data, formula, control = list())
```

## Arguments

- data:

  aggregated LFQData with a binary detection response

- formula:

  model formula with the binary response

- control:

  convergence controls, see
  [`strategy_logistf`](https://fgcz.github.io/prolfqua/reference/strategy.md)

## Value

a [ModelFirth](https://fgcz.github.io/prolfqua/reference/ModelFirth.md)
object holding one Firth logistic model per subject

## Note

Deprecated in prolfqua 1.9.5; use
`build_contrast_analysis(method = "firth")`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
