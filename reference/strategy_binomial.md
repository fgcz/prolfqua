# Create a quasibinomial detection-count strategy

Create a quasibinomial detection-count strategy

## Usage

``` r
strategy_binomial(
  modelstr,
  prior_count = 0.1,
  model_name = "binomial_nested",
  control = list()
)
```

## Arguments

- modelstr:

  right-hand-side model formula

- prior_count:

  non-negative symmetric pseudo-count

- model_name:

  model identity

- control:

  `epsilon` (default 1e-8) and `maxit` (default 25) of
  [`glm.control`](https://rdrr.io/r/stats/glm.control.html): the
  relative deviance change that stops IRLS and its iteration limit;
  relative, so independent of the number of detections, and glm's own
  defaults.

## Value

A
[`StrategyBinomial`](https://fgcz.github.io/prolfqua/reference/StrategyBinomial.md)
object.

## Examples

``` r
strategy <- strategy_binomial("~ group_", prior_count = 0.1)
strategy$formula
#> cbind(.detected, .undetected) ~ group_
#> <environment: 0x644bddc52668>
```
