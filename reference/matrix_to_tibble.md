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

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_uniprot_id_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`remove_na_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
x <- matrix(rnorm(20), ncol=4)
rownames(x) <- LETTERS[seq_len(nrow(x))]
matrix_to_tibble(x)
#> Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
#> `.name_repair` is omitted as of tibble 2.0.0.
#> ℹ Using compatibility `.name_repair`.
#> ℹ The deprecated feature was likely used in the prolfqua package.
#>   Please report the issue at <https://github.com/fgcz/prolfqua/issues>.
#> # A tibble: 5 × 5
#>   row.names       V1      V2      V3     V4
#>   <chr>        <dbl>   <dbl>   <dbl>  <dbl>
#> 1 A         -2.52    -1.61    0.545  -1.26 
#> 2 B         -1.49    -0.0719  1.06    0.707
#> 3 C         -1.38    -1.55    0.0597  0.508
#> 4 D         -0.511    0.0988 -0.359   1.62 
#> 5 E          0.00332  1.03   -0.764  -0.252
!(is.matrix(x) || is.data.frame(x))
#> [1] FALSE
```
