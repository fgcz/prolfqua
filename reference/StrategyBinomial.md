# Quasibinomial detection-count strategy

Quasibinomial detection-count strategy

Quasibinomial detection-count strategy

## Value

An R6 class generator.

## Details

Fits detected and undetected child-feature counts for each parent
feature using [`glm`](https://rdrr.io/r/stats/glm.html) with a
quasibinomial family. The symmetric pseudo-count stabilizes fits under
complete separation; it is not equivalent to Firth's bias-reducing
penalty.

## Public fields

- `formula`:

  quasibinomial model formula

- `model_name`:

  model identity

- `report_columns`:

  result columns supported by the strategy

- `is_mixed`:

  always FALSE

- `anova_df`:

  ANOVA extractor

- `prior_count`:

  symmetric pseudo-count added to both outcomes

## Methods

### Public methods

- [`StrategyBinomial$new()`](#method-StrategyBinomial-new)

- [`StrategyBinomial$model_fun()`](#method-StrategyBinomial-model_fun)

- [`StrategyBinomial$isSingular()`](#method-StrategyBinomial-isSingular)

- [`StrategyBinomial$contrast_fun()`](#method-StrategyBinomial-contrast_fun)

- [`StrategyBinomial$df_residual()`](#method-StrategyBinomial-df_residual)

- [`StrategyBinomial$sigma()`](#method-StrategyBinomial-sigma)

- [`StrategyBinomial$clone()`](#method-StrategyBinomial-clone)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

Create a quasibinomial count strategy.

#### Usage

    StrategyBinomial$new(
      modelstr,
      prior_count = 0.1,
      model_name = "binomial_nested",
      report_columns = c("statistic", "p.value", "p.value.adjusted", "moderated.p.value",
        "moderated.p.value.adjusted")
    )

#### Arguments

- `modelstr`:

  right-hand-side model formula, for example `"~ group_"`

- `prior_count`:

  non-negative symmetric pseudo-count

- `model_name`:

  model identity

- `report_columns`:

  result columns supported by the strategy

------------------------------------------------------------------------

### Method `model_fun()`

Fit one parent's detection counts.

#### Usage

    StrategyBinomial$model_fun(x, pb, get_formula = FALSE)

#### Arguments

- `x`:

  parent-by-sample count data

- `pb`:

  optional progress reporter

- `get_formula`:

  if TRUE, return the model formula without fitting

------------------------------------------------------------------------

### Method `isSingular()`

Check whether the model is singular.

#### Usage

    StrategyBinomial$isSingular(model)

#### Arguments

- `model`:

  fitted quasibinomial model

------------------------------------------------------------------------

### Method `contrast_fun()`

Compute linear contrasts.

#### Usage

    StrategyBinomial$contrast_fun(...)

#### Arguments

- `...`:

  passed to
  [`compute_contrast`](https://wolski.github.io/prolfqua/reference/compute_contrast.md)

------------------------------------------------------------------------

### Method `df_residual()`

Return residual degrees of freedom.

#### Usage

    StrategyBinomial$df_residual(model)

#### Arguments

- `model`:

  fitted quasibinomial model

------------------------------------------------------------------------

### Method [`sigma()`](https://rdrr.io/r/stats/sigma.html)

Return the Pearson residual scale used by
[`vcov()`](https://rdrr.io/r/stats/vcov.html).

#### Usage

    StrategyBinomial$sigma(model)

#### Arguments

- `model`:

  fitted quasibinomial model

------------------------------------------------------------------------

### Method `clone()`

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
