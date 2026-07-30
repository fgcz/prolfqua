# Create a quasibinomial detection-count strategy

Create a quasibinomial detection-count strategy

## Usage

``` r
strategy_binomial(
  modelstr,
  prior_count = 0.1,
  model_name = "binomial_nested",
  report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value",
    "moderated.p.value.adjusted")
)
```

## Arguments

- modelstr:

  right-hand-side model formula

- prior_count:

  non-negative symmetric pseudo-count

- model_name:

  model identity

- report_columns:

  result columns supported by the strategy

## Value

A
[`StrategyBinomial`](https://wolski.github.io/prolfqua/reference/StrategyBinomial.md)
object.

## Examples

``` r
strategy <- strategy_binomial("~ group_", prior_count = 0.1)
strategy$model_fun(get_formula = TRUE)
#> cbind(.detected, .undetected) ~ group_
#> <environment: 0x564609748ab8>
```
