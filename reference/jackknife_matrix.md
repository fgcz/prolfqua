# Compute correlation matrix with jack

Compute correlation matrix with jack

## Usage

``` r
jackknife_matrix(dataX, distmethod, ...)
```

## Arguments

- dataX:

  data.frame with transition intensities per peptide

- distmethod:

  dist or correlation method working with matrix i.e. cor

- ...:

  further parameters to method

## Value

summarizes results producced with jackknife

## See also

jackknife

Other transitioncorrlation:
[`cor_jackknife_matrix()`](https://wolski.github.io/prolfqua/reference/cor_jackknife_matrix.md),
[`cor_order()`](https://wolski.github.io/prolfqua/reference/cor_order.md),
[`jackknife()`](https://wolski.github.io/prolfqua/reference/jackknife.md)

## Examples

``` r
dataX <- matrix(rnorm(20), ncol=4)
rownames(dataX)<- paste("R",seq_len(nrow(dataX)),sep="")
colnames(dataX)<- paste("C",seq_len(ncol(dataX)),sep="")
tmp <- jackknife(dataX, cor, use="pairwise.complete.obs", method="pearson")
res <- jackknife_matrix(dataX, cor)
res
#>           C1         C2         C3          C4
#> C1 1.0000000  0.7068741 0.73327820  0.58869599
#> C2 0.7068741  1.0000000 0.28845670 -0.17789058
#> C3 0.7332782  0.2884567 1.00000000  0.04364043
#> C4 0.5886960 -0.1778906 0.04364043  1.00000000
stopifnot(dim(res) == c(4,4))
res <- jackknife_matrix(dataX, cor, method="spearman")
stopifnot(dim(res) == c(4,4))
```
