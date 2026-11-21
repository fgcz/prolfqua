# Create a Firth-penalized quasibinomial detection-count strategy

Create a Firth-penalized quasibinomial detection-count strategy

## Usage

``` r
strategy_binomial_firth(
  modelstr,
  model_name = "binomial_firth_nested",
  control = list()
)
```

## Arguments

- modelstr:

  right-hand-side model formula

- model_name:

  model identity

- control:

  `maxit` (default 100) and `tol` (default 1e-8): iteration limit of the
  Firth adjustment and the largest coefficient change (log odds) that
  stops it; log odds of detection range over about -10 to 10, so 1e-8 is
  far below any reported precision.

## Value

A
[`StrategyBinomialFirth`](https://fgcz.github.io/prolfqua/reference/StrategyBinomialFirth.md)
object.

## Examples

``` r
strategy <- strategy_binomial_firth("~ group_")
strategy$formula
#> cbind(.detected, .undetected) ~ group_
#> <environment: 0x59a8db211cc0>
```
