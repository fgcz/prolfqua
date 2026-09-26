# matrix or data.frame to tibble (taken from tidyquant)

matrix or data.frame to tibble (taken from tidyquant)

## Usage

``` r
matrix_to_tibble(x, preserve_row_names = "row.names", ...)
```

## Arguments

- x:

  a matrix

- preserve_row_names:

  give name to rownames column, if NULL discard rownames

- ...:

  further parameters passed to as_tibble

## Value

The computed result.

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_uniprot_id_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`is_contaminant()`](https://wolski.github.io/prolfqua/reference/is_contaminant.md),
[`is_decoy()`](https://wolski.github.io/prolfqua/reference/is_decoy.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
x <- matrix(rnorm(20), ncol=4)
rownames(x) <- LETTERS[seq_len(nrow(x))]
matrix_to_tibble(x)
#> # A tibble: 5 × 5
#>   row.names     V1    V2     V3     V4
#>   <chr>      <dbl> <dbl>  <dbl>  <dbl>
#> 1 A          0.344 0.841 -1.05   1.76 
#> 2 B         -0.830 1.45   1.84   0.476
#> 3 C          0.482 0.137 -2.23   0.783
#> 4 D         -0.689 0.326 -0.596 -1.15 
#> 5 E         -1.09  1.22  -0.230  0.950
!(is.matrix(x) || is.data.frame(x))
#> [1] FALSE
```
