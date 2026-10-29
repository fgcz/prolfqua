# Create a Firth-penalized quasibinomial detection-count strategy

Create a Firth-penalized quasibinomial detection-count strategy

## Usage

``` r
strategy_binomial_firth(modelstr, model_name = "binomial_firth_nested")
```

## Arguments

- modelstr:

  right-hand-side model formula

- model_name:

  model identity

## Value

A
[`StrategyBinomialFirth`](https://fgcz.github.io/prolfqua/reference/StrategyBinomialFirth.md)
object.

## Examples

``` r
strategy <- strategy_binomial_firth("~ group_")
strategy$formula
#> cbind(.detected, .undetected) ~ group_
#> <environment: 0x5645ed68bf50>
```
