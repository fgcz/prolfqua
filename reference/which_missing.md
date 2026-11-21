# Draw missing values for an intensity vector based on its values

Draw missing values for an intensity vector based on its values

## Usage

``` r
which_missing(x, weight_missing = 0.2)
```

## Arguments

- x:

  vector of intensities

- weight_missing:

  greater weight more missing

## Value

A logical vector indicating which values are missing.

## Examples

``` r
which_missing(2**rnorm(10,2,0.4))
#>  [1] FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE
```
