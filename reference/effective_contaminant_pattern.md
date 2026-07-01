# Effective contaminant regex actually applied by [`is_contaminant`](https://wolski.github.io/prolfqua/reference/is_contaminant.md)

Effective contaminant regex actually applied by
[`is_contaminant`](https://wolski.github.io/prolfqua/reference/is_contaminant.md)

## Usage

``` r
effective_contaminant_pattern(pattern = NULL)
```

## Arguments

- pattern:

  optional configured contaminant regex

## Value

a single regex string

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`effective_decoy_pattern()`](https://wolski.github.io/prolfqua/reference/effective_decoy_pattern.md),
[`get_uniprot_id_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`is_contaminant()`](https://wolski.github.io/prolfqua/reference/is_contaminant.md),
[`is_decoy()`](https://wolski.github.io/prolfqua/reference/is_decoy.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`remove_na_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
effective_contaminant_pattern()
#> [1] "^zz|^CON__|^CON_|^Cont_|^contam_"
```
