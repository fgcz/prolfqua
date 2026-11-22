# R6 interface class representing modelling result

R6 interface class representing modelling result

## Value

An R6 class generator.

## Methods

### Public methods

- [`ModelInterface$provenance()`](#method-ModelInterface-provenance)

- [`ModelInterface$get_provenance_log()`](#method-ModelInterface-get_provenance_log)

- [`ModelInterface$get_coefficients()`](#method-ModelInterface-get_coefficients)

- [`ModelInterface$get_anova()`](#method-ModelInterface-get_anova)

- [`ModelInterface$coef_histogram()`](#method-ModelInterface-coef_histogram)

- [`ModelInterface$coef_volcano()`](#method-ModelInterface-coef_volcano)

- [`ModelInterface$coef_pairs()`](#method-ModelInterface-coef_pairs)

- [`ModelInterface$anova_histogram()`](#method-ModelInterface-anova_histogram)

- [`ModelInterface$clone()`](#method-ModelInterface-clone)

------------------------------------------------------------------------

### `ModelInterface$provenance()`

problems recorded while fitting: failed fits, fit warnings (e.g.
non-convergence), singular fits, fits refit at the limit of detection,
missing responses dropped

#### Usage

    ModelInterface$provenance()

#### Returns

tibble with columns step, event, n, ids (list), reason, failure,
response

------------------------------------------------------------------------

### `ModelInterface$get_provenance_log()`

internal: the model's provenance record

#### Usage

    ModelInterface$get_provenance_log()

#### Returns

Provenance object

------------------------------------------------------------------------

### `ModelInterface$get_coefficients()`

return model coefficients

#### Usage

    ModelInterface$get_coefficients()

#### Returns

data.frame

------------------------------------------------------------------------

### `ModelInterface$get_anova()`

perform ANOVA analysis

#### Usage

    ModelInterface$get_anova()

#### Returns

data.frame

------------------------------------------------------------------------

### `ModelInterface$coef_histogram()`

histogram of coefficient p-values

#### Usage

    ModelInterface$coef_histogram()

#### Returns

list with `plot` (ggplot) and file `name`

------------------------------------------------------------------------

### `ModelInterface$coef_volcano()`

volcano plot of non-intercept coefficients

#### Usage

    ModelInterface$coef_volcano()

#### Returns

list with `plot` (ggplot) and file `name`

------------------------------------------------------------------------

### `ModelInterface$coef_pairs()`

coefficient estimates in wide format, one column per coefficient, for a
pairs plot

#### Usage

    ModelInterface$coef_pairs()

#### Returns

list with `plot` (data.frame) and file `name`

------------------------------------------------------------------------

### `ModelInterface$anova_histogram()`

histogram of ANOVA p-values or FDR

#### Usage

    ModelInterface$anova_histogram(what = c("p.value", "FDR"))

#### Arguments

- `what`:

  show either "p.value" or "FDR"

#### Returns

list with `plot` (ggplot) and file `name`

------------------------------------------------------------------------

### `ModelInterface$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ModelInterface$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
mi <- ModelInterface$new()

abstract <- c("get_coefficients", "get_anova", "coef_histogram", "coef_volcano", "coef_pairs",
  "anova_histogram")
for (m in abstract) {
  stopifnot(inherits(try(mi[[m]](), silent = TRUE), "try-error"))
}

```
