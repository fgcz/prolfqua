# add missing values to x vector based on the values of x

add missing values to x vector based on the values of x

## Usage

``` r
which_missing(x, weight_missing = 0.2)
```

## Arguments

- x:

  vector of intensities

- weight_missing:

  greater weight more missing

## Examples

``` r
which_missing(2**rnorm(10,2,0.4))
#>  [1] FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE
```
