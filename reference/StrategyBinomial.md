# Quasibinomial detection-count strategy

Fits detected and undetected child-feature counts for each parent
feature using [`glm`](https://rdrr.io/r/stats/glm.html) with a
quasibinomial family. The symmetric pseudo-count stabilizes fits under
complete separation; it is not equivalent to Firth's bias-reducing
penalty.

## Value

An R6 class generator.

## Super class

[`StrategyBase`](https://fgcz.github.io/prolfqua/reference/StrategyBase.md)
-\> `StrategyBinomial`

## Public fields

- `prior_count`:

  symmetric pseudo-count added to both outcomes

## Methods

### Public methods

- [`StrategyBinomial$new()`](#method-StrategyBinomial-initialize)

- [`StrategyBinomial$clone()`](#method-StrategyBinomial-clone)

Inherited methods

- [`StrategyBase$contrast_fun()`](https://fgcz.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-contrast_fun)
- [`StrategyBase$df_residual()`](https://fgcz.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-df_residual)
- [`StrategyBase$isSingular()`](https://fgcz.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-isSingular)
- [`StrategyBase$model_fun()`](https://fgcz.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-model_fun)
- [`StrategyBase$sigma()`](https://fgcz.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-sigma)

------------------------------------------------------------------------

### `StrategyBinomial$new()`

Create a quasibinomial count strategy.
[`sigma()`](https://rdrr.io/r/stats/sigma.html) is the Pearson residual
scale used by [`vcov()`](https://rdrr.io/r/stats/vcov.html).

#### Usage

    StrategyBinomial$new(
      modelstr,
      prior_count = 0.1,
      model_name = "binomial_nested"
    )

#### Arguments

- `modelstr`:

  right-hand-side model formula, for example `"~ group_"`

- `prior_count`:

  non-negative symmetric pseudo-count

- `model_name`:

  model identity

------------------------------------------------------------------------

### `StrategyBinomial$clone()`

The objects of this class are cloneable with this method.

#### Usage

    StrategyBinomial$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
dat <- data.frame(
  group_ = factor(rep(c("A", "B"), each = 4)),
  detected = c(1, 2, 1, 3, 4, 5, 3, 5),
  undetected = c(4, 3, 4, 2, 1, 0, 2, 0)
)
strategy <- StrategyBinomial$new("~ group_")
fit <- strategy$model_fun(dat)
coefficients(fit)
#> (Intercept)     group_B 
#>  -0.5937747   2.2264695 
```
