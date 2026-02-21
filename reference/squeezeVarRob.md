# Robustly Squeeze Sample Variances

This function squeezes a set of sample variances together by computing
empirical Bayes posterior means in a way that is robust against the
presence of very small non-integer degrees of freedom values.

## Usage

``` r
squeezeVarRob(
  var,
  df,
  covariate = NULL,
  robust = FALSE,
  winsor.tail.p = c(0.05, 0.1),
  min_df = 1
)
```

## Arguments

- var:

  A numeric vector of independent sample variances.

- df:

  A numeric vector of degrees of freedom for the sample variances.

- covariate:

  If `non-NULL`, `var.prior` will depend on this numeric covariate.
  Otherwise, `var.prior` is constant.

- robust:

  A logical indicating wheter the estimation of `df.prior` and
  `var.prior` should be robustified against outlier sample variances.
  Defaults to `FALSE`.

- winsor.tail.p:

  A numeric vector of length 1 or 2, giving left and right tail
  proportions of `x` to Winsorize. Used only when `robust=TRUE`.

- min_df:

  A numeric value indicating the minimal degrees of freedom that will be
  taken into account for calculating the prior degrees of freedom and
  prior variance.

## Value

A list with components: `var.post` A numeric vector of posterior
variances. `var.prior` The location of prior distribution. A vector if
`covariate` is non-`NULL`, otherwise a scalar. `df.prior` The degrees of
freedom of prior distribution. A vector if `robust=TRUE`, otherwise a
scalar.

## Examples

``` r
var <- rexp(1000)
df <- sample( 3:10, 1000, replace = TRUE)
tmp <- squeezeVarRob(var, df)
tmp <- squeezeVarRob(var, df, robust = TRUE)
```
