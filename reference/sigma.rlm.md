# Residual scale estimate for rlm objects

[`stats::sigma`](https://rdrr.io/r/stats/sigma.html) returns `NA` for
[`rlm`](https://rdrr.io/pkg/MASS/man/rlm.html) objects. This S3 method
computes a weighted residual scale estimate instead.

## Usage

``` r
# S3 method for class 'rlm'
sigma(object, ...)
```

## Arguments

- object:

  an `rlm` object

- ...:

  ignored

## Value

numeric scalar
