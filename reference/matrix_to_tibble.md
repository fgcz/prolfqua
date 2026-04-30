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
#>   row.names       V1      V2      V3      V4
#>   <chr>        <dbl>   <dbl>   <dbl>   <dbl>
#> 1 A          0.0324  -1.81   -0.103   0.972 
#> 2 B         -0.00624  1.25   -0.0252  0.0103
#> 3 C         -0.505   -0.651   1.01   -0.598 
#> 4 D          0.471   -0.0453  0.223   0.621 
#> 5 E          1.56     1.15    0.747  -0.107 
!(is.matrix(x) || is.data.frame(x))
#> [1] FALSE
```
