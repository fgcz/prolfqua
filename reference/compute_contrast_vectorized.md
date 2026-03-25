# Vectorized version of [`compute_contrast`](https://wolski.github.io/prolfqua/reference/compute_contrast.md)

Same semantics but uses vectorized matrix multiplication instead of a
per-row loop. NAs in `coefficients(m)` propagate naturally via
`linfct %*% coef` for rows that reference missing coefficients.

## Usage

``` r
compute_contrast_vectorized(m, linfct, confint = 0.95)
```

## Arguments

- m:

  linear model generated using lm

- linfct:

  linear function

- confint:

  confidence interval default 0.95
