# R6 interface class representing modelling result

R6 interface class representing modelling result

R6 interface class representing modelling result

## Value

An R6 class generator.

## Methods

### Public methods

- [`ModelInterface$get_coefficients()`](#method-ModelInterface-get_coefficients)

- [`ModelInterface$get_anova()`](#method-ModelInterface-get_anova)

- [`ModelInterface$coef_histogram()`](#method-ModelInterface-coef_histogram)

- [`ModelInterface$coef_volcano()`](#method-ModelInterface-coef_volcano)

- [`ModelInterface$coef_pairs()`](#method-ModelInterface-coef_pairs)

- [`ModelInterface$anova_histogram()`](#method-ModelInterface-anova_histogram)

- [`ModelInterface$clone()`](#method-ModelInterface-clone)

------------------------------------------------------------------------

### Method `get_coefficients()`

return model coefficients

#### Usage

    ModelInterface$get_coefficients()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `get_anova()`

perform ANOVA analysis

#### Usage

    ModelInterface$get_anova()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `coef_histogram()`

histogram of coefficient p-values

#### Usage

    ModelInterface$coef_histogram()

#### Returns

list with `plot` (ggplot) and file `name`

------------------------------------------------------------------------

### Method `coef_volcano()`

volcano plot of non-intercept coefficients

#### Usage

    ModelInterface$coef_volcano()

#### Returns

list with `plot` (ggplot) and file `name`

------------------------------------------------------------------------

### Method `coef_pairs()`

coefficient estimates in wide format, one column per coefficient, for a
pairs plot

#### Usage

    ModelInterface$coef_pairs()

#### Returns

list with `plot` (data.frame) and file `name`

------------------------------------------------------------------------

### Method `anova_histogram()`

histogram of ANOVA p-values or FDR

#### Usage

    ModelInterface$anova_histogram(what = c("p.value", "FDR"))

#### Arguments

- `what`:

  show either "p.value" or "FDR"

#### Returns

list with `plot` (ggplot) and file `name`

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ModelInterface$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
mi <- ModelInterface$new()

testthat::expect_error(mi$get_coefficients())
testthat::expect_error(mi$get_anova())
testthat::expect_error(mi$coef_histogram())
testthat::expect_error(mi$coef_volcano())
testthat::expect_error(mi$coef_pairs())
testthat::expect_error(mi$anova_histogram())

```
