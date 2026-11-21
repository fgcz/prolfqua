# Correlation panel for pairs plots (default lower panel of `pairs_smooth()`)

Correlation panel for pairs plots (default lower panel of
[`pairs_smooth()`](https://fgcz.github.io/prolfqua/reference/pairs_smooth.md))

## Usage

``` r
panel_cor(x, y, digits = 2, ...)
```

## Arguments

- x:

  numeric data

- y:

  numeric data

- digits:

  number of digits to display

- ...:

  not used

## Value

`NULL`, invisibly; called for its side effect of writing the
correlation, R^2 and p-value into the current panel

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_uniprot_id_from_fasta_header()`](https://fgcz.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`is_contaminant()`](https://fgcz.github.io/prolfqua/reference/is_contaminant.md),
[`is_decoy()`](https://fgcz.github.io/prolfqua/reference/is_decoy.md),
[`matrix_to_tibble()`](https://fgcz.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://fgcz.github.io/prolfqua/reference/multigroup_volcano.md),
[`pairs_smooth()`](https://fgcz.github.io/prolfqua/reference/pairs_smooth.md),
[`table_facade()`](https://fgcz.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
graphics::plot(1:3, 1:3, type = "n")
panel_cor(1:3, 1:3)
```
