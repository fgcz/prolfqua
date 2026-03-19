# R6 interface class representing modelling result

R6 interface class representing modelling result

R6 interface class representing modelling result

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

plot histogram of coefficients

#### Usage

    ModelInterface$coef_histogram()

#### Returns

ggplot

------------------------------------------------------------------------

### Method `coef_volcano()`

plot volcano of coefficients

#### Usage

    ModelInterface$coef_volcano()

#### Returns

ggplot

------------------------------------------------------------------------

### Method `coef_pairs()`

pairs plot of coefficients

#### Usage

    ModelInterface$coef_pairs()

#### Returns

ggplot

------------------------------------------------------------------------

### Method `anova_histogram()`

histogram of p-values and FDR for anova results

#### Usage

    ModelInterface$anova_histogram()

#### Returns

ggplot

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
