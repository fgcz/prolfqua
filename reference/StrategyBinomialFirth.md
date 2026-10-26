# Firth-penalized quasibinomial detection-count strategy

Firth-penalized quasibinomial detection-count strategy

Firth-penalized quasibinomial detection-count strategy

## Value

An R6 class generator.

## Details

Fits detected and undetected child-feature counts like
[`StrategyBinomial`](https://wolski.github.io/prolfqua/reference/StrategyBinomial.md),
but replaces the fixed pseudo-count with Firth's bias reduction: each
count row gets \\h/2\\ added to both outcomes, where \\h\\ is its hat
value at the fitted probabilities, iterated to convergence. The
coefficients are Firth's penalized maximum-likelihood estimates, finite
under complete separation. [`vcov()`](https://rdrr.io/r/stats/vcov.html)
and [`sigma()`](https://rdrr.io/r/stats/sigma.html) are those of the
quasibinomial fit to the adjusted counts, so the dispersion is estimated
from the scatter of the counts between samples.

## Super classes

[`prolfqua::StrategyBase`](https://wolski.github.io/prolfqua/reference/StrategyBase.md)
-\>
[`prolfqua::StrategyBinomial`](https://wolski.github.io/prolfqua/reference/StrategyBinomial.md)
-\> `StrategyBinomialFirth`

## Methods

### Public methods

- [`StrategyBinomialFirth$new()`](#method-StrategyBinomialFirth-new)

- [`StrategyBinomialFirth$clone()`](#method-StrategyBinomialFirth-clone)

Inherited methods

- [`prolfqua::StrategyBase$contrast_fun()`](https://wolski.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-contrast_fun)
- [`prolfqua::StrategyBase$df_residual()`](https://wolski.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-df_residual)
- [`prolfqua::StrategyBase$isSingular()`](https://wolski.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-isSingular)
- [`prolfqua::StrategyBase$model_fun()`](https://wolski.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-model_fun)
- [`prolfqua::StrategyBase$sigma()`](https://wolski.github.io/prolfqua/html/StrategyBase.html#method-StrategyBase-sigma)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

Create a Firth-penalized quasibinomial count strategy.

#### Usage

    StrategyBinomialFirth$new(modelstr, model_name = "binomial_firth_nested")

#### Arguments

- `modelstr`:

  right-hand-side model formula, for example `"~ group_"`

- `model_name`:

  model identity

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    StrategyBinomialFirth$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
dat <- data.frame(
  group_ = factor(rep(c("A", "B"), each = 4)),
  detected = c(0, 0, 0, 0, 4, 5, 3, 5),
  undetected = c(5, 5, 5, 5, 1, 0, 2, 0)
)
strategy <- StrategyBinomialFirth$new("~ group_")
fit <- strategy$model_fun(dat)
coefficients(fit)
#> (Intercept)     group_B 
#>   -3.713572    5.323010 
```
