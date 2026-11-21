# Robustly scale a matrix by medians and MADs

Robustly scale a matrix by medians and MADs

## Usage

``` r
robust_scale(data, dim = 2, preserve_mean = FALSE)
```

## Value

a matrix like `data`, centred by the median of each column (`dim = 2`)
or row (`dim = 1`) and divided by its MAD relative to the mean MAD; only
centred, with a warning, if a MAD is zero; with `preserve_mean = TRUE`
the mean of the medians is added back

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://fgcz.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://fgcz.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://fgcz.github.io/prolfqua/reference/get_robscales.md),
[`nr_B_in_A()`](https://fgcz.github.io/prolfqua/reference/nr_B_in_A.md),
[`scale_with_subset()`](https://fgcz.github.io/prolfqua/reference/scale_with_subset.md)

## Examples

``` r
mat <- matrix(c(1, 2, 3, 4, 10, 12), nrow = 3)
scaled <- robust_scale(mat)
dim(scaled)
#> [1] 3 2
```
