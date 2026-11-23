# Draw missing values for an intensity vector based on its values

Each value goes missing at random with a probability that falls with its
intensity (through its normal quantile), so low intensities go missing
more often, as in LFQ data.

## Usage

``` r
which_missing(x, weight_missing = 0.2)
```

## Arguments

- x:

  numeric vector of intensities without NA

- weight_missing:

  one number in \[0, 1\]; greater weight means more missing values

## Value

A logical vector indicating which values are missing.

## Examples

``` r
which_missing(2**rnorm(10,2,0.4))
#>  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
```
