# Base class for per-subject model strategies

Holds the fields, constructor and default methods shared by
[`StrategyLM`](https://fgcz.github.io/prolfqua/reference/StrategyLM.md),
[`StrategyRLM`](https://fgcz.github.io/prolfqua/reference/StrategyRLM.md),
[`StrategyLmer`](https://fgcz.github.io/prolfqua/reference/StrategyLmer.md),
[`StrategyRfit`](https://fgcz.github.io/prolfqua/reference/StrategyRfit.md),
[`StrategyLogistf`](https://fgcz.github.io/prolfqua/reference/StrategyLogistf.md)
and
[`StrategyBinomial`](https://fgcz.github.io/prolfqua/reference/StrategyBinomial.md).
A subclass implements the private `fit(x)` method and overrides only the
methods that differ.

## Value

An R6 class generator.

## Public fields

- `formula`:

  model formula

- `model_name`:

  name of model

- `anova_df`:

  ANOVA extractor, see
  [`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md)

## Methods

### Public methods

- [`StrategyBase$new()`](#method-StrategyBase-initialize)

- [`StrategyBase$model_fun()`](#method-StrategyBase-model_fun)

- [`StrategyBase$isSingular()`](#method-StrategyBase-isSingular)

- [`StrategyBase$contrast_fun()`](#method-StrategyBase-contrast_fun)

- [`StrategyBase$df_residual()`](#method-StrategyBase-df_residual)

- [`StrategyBase$sigma()`](#method-StrategyBase-sigma)

- [`StrategyBase$clone()`](#method-StrategyBase-clone)

------------------------------------------------------------------------

### `StrategyBase$new()`

Create a new strategy

#### Usage

    StrategyBase$new(modelstr, model_name = "Model")

#### Arguments

- `modelstr`:

  model formula string

- `model_name`:

  name of model

------------------------------------------------------------------------

### `StrategyBase$model_fun()`

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

### `StrategyBase$isSingular()`

Check if model is singular (NA coefficients or df \< 2)

#### Usage

    StrategyBase$isSingular(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### `StrategyBase$contrast_fun()`

Compute contrasts from fitted model

#### Usage

    StrategyBase$contrast_fun(...)

#### Arguments

- `...`:

  passed to
  [`compute_contrast`](https://fgcz.github.io/prolfqua/reference/compute_contrast.md)

------------------------------------------------------------------------

### `StrategyBase$df_residual()`

Get residual degrees of freedom

#### Usage

    StrategyBase$df_residual(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### `StrategyBase$sigma()`

Get residual standard error

#### Usage

    StrategyBase$sigma(model)

#### Arguments

- `model`:

  fitted model

------------------------------------------------------------------------

### `StrategyBase$clone()`

The objects of this class are cloneable with this method.

#### Usage

    StrategyBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
