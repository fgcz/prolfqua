# Get group-level linear functions from a fitted model

Get group-level linear functions from a fitted model

## Usage

``` r
linfct_from_model(m, as_list = TRUE)
```

## Arguments

- m:

  linear model

- as_list:

  if TRUE, return the factor and interaction matrices as a list; if
  FALSE, as one matrix

## Value

if `as_list = TRUE`, a list with the matrices `linfct_factors` (one row
per factor level) and `linfct_interactions` (one row per observed factor
combination), each with one column per model coefficient; otherwise the
two matrices bound by rows into a single matrix

## Details

The fitted model supplies its coefficient names and observed categorical
factor combinations. Each combination is expressed as weights on the
fitted coefficients; weights for an individual factor level are averaged
equally over the observed combinations containing that level. Thus
`linfct_factors` describes marginal factor-level means and
`linfct_interactions` describes cell means. With treatment coding, a
reference-level mean includes the intercept; a non-reference mean also
includes its treatment coefficient. The *Under the hood* article shows a
complete `~ A * B` example. These helpers assume categorical predictors
and the usual intercept-based coding; they are not a general
marginal-means implementation for continuous predictors or arbitrary
coding.

## Note

Deprecated in prolfqua 1.9.5; use `Contrasts$get_linfct()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
