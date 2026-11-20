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

- [`StrategyBinomial$zero_residual_variance()`](#method-StrategyBinomial-zero_residual_variance)

- [`StrategyBinomial$clone()`](#method-StrategyBinomial-clone)

Inherited methods

- [`StrategyBase$contrast_fun()`](https://fgcz.github.io/prolfqua/reference/StrategyBase.html#method-contrast_fun)
- [`StrategyBase$df_residual()`](https://fgcz.github.io/prolfqua/reference/StrategyBase.html#method-df_residual)
- [`StrategyBase$isSingular()`](https://fgcz.github.io/prolfqua/reference/StrategyBase.html#method-isSingular)
- [`StrategyBase$model_fun()`](https://fgcz.github.io/prolfqua/reference/StrategyBase.html#method-model_fun)
- [`StrategyBase$sigma()`](https://fgcz.github.io/prolfqua/reference/StrategyBase.html#method-sigma)

------------------------------------------------------------------------

### `StrategyBinomial$new()`

Create a quasibinomial count strategy.
[`sigma()`](https://rdrr.io/r/stats/sigma.html) is the Pearson residual
scale used by [`vcov()`](https://rdrr.io/r/stats/vcov.html).

#### Usage

    StrategyBinomial$new(
      modelstr,
      prior_count = 0.1,
      model_name = "binomial_nested",
      control = list()
    )

#### Arguments

- `modelstr`:

  right-hand-side model formula, for example `"~ group_"`

- `prior_count`:

  non-negative symmetric pseudo-count

- `model_name`:

  model identity

- `control`:

  named list of convergence controls, see
  [`strategy_binomial`](https://fgcz.github.io/prolfqua/reference/strategy_binomial.md)

------------------------------------------------------------------------

### `StrategyBinomial$zero_residual_variance()`

Never flags: a zero Pearson dispersion (identical detection pattern
within groups) is bounded below by the binomial variance
(`binomial_bound`), which keeps the test defined

#### Usage

    StrategyBinomial$zero_residual_variance(sigma, scale)

#### Arguments

- `sigma`:

  Pearson residual scales

- `scale`:

  unused

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
