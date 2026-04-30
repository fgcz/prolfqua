# robust scale wrapper

robust scale wrapper

## Usage

``` r
robust_scale(data, dim = 2, preserve_mean = FALSE)
```

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://wolski.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md),
[`get_robscales()`](https://wolski.github.io/prolfqua/reference/get_robscales.md),
[`scale_with_subset()`](https://wolski.github.io/prolfqua/reference/scale_with_subset.md)

## Examples

``` r
mat <- matrix(c(1, 2, 3, 4, 10, 12), nrow = 3)
scaled <- robust_scale(mat)
dim(scaled)
#> [1] 3 2
```
