# Create a quasibinomial detection-count strategy

Create a quasibinomial detection-count strategy

## Usage

``` r
strategy_binomial(modelstr, prior_count = 0.1, model_name = "binomial_nested")
```

## Arguments

- modelstr:

  right-hand-side model formula

- prior_count:

  non-negative symmetric pseudo-count

- model_name:

  model identity

## Value

A
[`StrategyBinomial`](https://wolski.github.io/prolfqua/reference/StrategyBinomial.md)
object.

## Examples

``` r
strategy <- strategy_binomial("~ group_", prior_count = 0.1)
strategy$formula
#> cbind(.detected, .undetected) ~ group_
#> <environment: 0x5593e90387b8>
```
