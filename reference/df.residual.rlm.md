# Residual degrees of freedom for rlm objects

[`stats::df.residual`](https://rdrr.io/r/stats/df.residual.html) returns
`NA` for [`rlm`](https://rdrr.io/pkg/MASS/man/rlm.html) objects. This S3
method computes weighted residual df instead.

## Usage

``` r
# S3 method for class 'rlm'
df.residual(object, ...)
```

## Arguments

- object:

  an `rlm` object

- ...:

  ignored

## Value

numeric scalar
