# Compute jackknive resampling on matrix

Compute jackknive resampling on matrix

## Usage

``` r
jackknife(xdata, .method, ...)
```

## Arguments

- xdata:

  matrix

- .method:

  method i.e. cor, parameters

- ...:

  further parameters to .method

## Value

list with all jackknife matrices

list of matrices

## See also

Other transitioncorrlation:
[`cor_jackknife_matrix()`](https://wolski.github.io/prolfqua/reference/cor_jackknife_matrix.md),
[`cor_order()`](https://wolski.github.io/prolfqua/reference/cor_order.md),
[`jackknife_matrix()`](https://wolski.github.io/prolfqua/reference/jackknife_matrix.md)

## Examples

``` r
xx <- matrix(rnorm(20), ncol=4)
cortest <- function(x){print(dim(x));cor(x)}
jackknife(xx, cortest)
#> [1] 4 4
#> [1] 4 4
#> [1] 4 4
#> [1] 4 4
#> [1] 4 4
#> [1] 5 4
tmp <- jackknife(xx, cor, use="pairwise.complete.obs", method="pearson")
stopifnot(length(tmp) == 3)
stopifnot(length(tmp[[2]]) == nrow(xx))
```
