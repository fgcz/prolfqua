# Median polish estimate of e.g. protein from peptide intensities

Compute Tukey's median polish estimate of a protein from peptide or
precursor intensities

## Usage

``` r
medpolish_estimate(x, name = FALSE, sampleName = "sampleName")
```

## Arguments

- x:

  a matrix

- name:

  if TRUE returns the name of the summary column

## Value

data.frame with number of rows equal to number of columns of input
matrix.

## See also

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_topN()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_topN.md),
[`estimate_intensity()`](https://wolski.github.io/prolfqua/reference/estimate_intensity.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`medpolish_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`rlm_estimate()`](https://wolski.github.io/prolfqua/reference/rlm_estimate.md),
[`rlm_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/rlm_estimate_dfconfig.md)

## Examples

``` r
medpolish_estimate(name = TRUE)
#> [1] "medpolish"
gg <- matrix(runif(20), 4, 5)
rownames(gg) <- paste0("A", 1:4)
colnames(gg) <- make.names(1:5)
gg
#>            X1        X2        X3        X4        X5
#> A1 0.65767434 0.1729643 0.2074822 0.2092584 0.5467840
#> A2 0.48037155 0.4603666 0.1351013 0.7143354 0.7838502
#> A3 0.01651635 0.5609443 0.2820738 0.5724989 0.7168003
#> A4 0.96492654 0.7814906 0.4122952 0.3613100 0.7561737
mx <- medpolish_estimate(gg)
```
