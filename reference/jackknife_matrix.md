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
#>            C1         C2         C3         C4
#> C1  1.0000000 -0.1838474 -0.5071304  0.9903309
#> C2 -0.1838474  1.0000000  0.7036950 -0.7823525
#> C3 -0.5071304  0.7036950  1.0000000  0.3141103
#> C4  0.9903309 -0.7823525  0.3141103  1.0000000
stopifnot(dim(res) == c(4,4))
res <- jackknife_matrix(dataX, cor, method="spearman")
stopifnot(dim(res) == c(4,4))
```
