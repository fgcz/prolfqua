# Base class for per-subject model strategies

Base class for per-subject model strategies

Base class for per-subject model strategies

## Value

An R6 class generator.

## Details

Holds the fields, constructor and default methods shared by
[`StrategyLM`](https://wolski.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyRLM`](https://wolski.github.io/prolfqua/reference/StrategyRLM.md),
[`StrategyLmer`](https://wolski.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyRfit`](https://wolski.github.io/prolfqua/reference/StrategyRfit.md),
[`StrategyLogistf`](https://wolski.github.io/prolfqua/reference/StrategyLogistf.md)
and
[`StrategyBinomial`](https://wolski.github.io/prolfqua/reference/StrategyBinomial.md).
A subclass implements the private `fit(x)` method and overrides only the
methods that differ.

## Public fields

- `formula`:

  model formula

- `model_name`:

  name of model

- `anova_df`:

  ANOVA extractor, see
  [`AnovaExtractor`](https://wolski.github.io/prolfqua/reference/AnovaExtractor.md)

## Methods

### Public methods

- [`StrategyBase$new()`](#method-StrategyBase-new)

- [`StrategyBase$model_fun()`](#method-StrategyBase-model_fun)

- [`StrategyBase$isSingular()`](#method-StrategyBase-isSingular)

- [`StrategyBase$contrast_fun()`](#method-StrategyBase-contrast_fun)

- [`StrategyBase$df_residual()`](#method-StrategyBase-df_residual)

- [`StrategyBase$sigma()`](#method-StrategyBase-sigma)

- [`StrategyBase$clone()`](#method-StrategyBase-clone)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

Create a new strategy

#### Usage

    StrategyBase$new(modelstr, model_name = "Model")

#### Arguments

- `modelstr`:

  model formula string

- `model_name`:

  name of model

------------------------------------------------------------------------

### Method `model_fun()`

Fit the model to one subject's data. A failed fit returns the error
message.

#### Usage

    StrategyBase$model_fun(x, pb)

#### Arguments

- `x`:

  data.frame for one subject

- `pb`:

  optional progress bar

------------------------------------------------------------------------

### Method `isSingular()`

Check if model is singular (NA coefficients or df \< 2)

#### Usage

    StrategyBase$isSingular(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### Method `contrast_fun()`

Compute contrasts from fitted model

#### Usage

    StrategyBase$contrast_fun(...)

#### Arguments

- `...`:

  passed to
  [`compute_contrast`](https://wolski.github.io/prolfqua/reference/compute_contrast.md)

------------------------------------------------------------------------

### Method `df_residual()`

Get residual degrees of freedom

#### Usage

    StrategyBase$df_residual(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### Method [`sigma()`](https://rdrr.io/r/stats/sigma.html)

Get residual standard error

#### Usage

    StrategyBase$sigma(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    StrategyBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
