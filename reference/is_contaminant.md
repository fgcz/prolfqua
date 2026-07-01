# Detect contaminant identifiers

As
[`is_decoy`](https://wolski.github.io/prolfqua/reference/is_decoy.md),
but for contaminant entries (keratin, trypsin, BSA, ...). Built-in
anchored default prefixes unioned with an optional configured `pattern`;
empty / `NULL` uses the defaults only.

## Usage

``` r
is_contaminant(ids, pattern = NULL)
```

## Arguments

- ids:

  character vector of (prefixed) identifiers

- pattern:

  optional additional contaminant regex, unioned with defaults

## Value

logical vector, `TRUE` where an id looks like a contaminant

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`effective_contaminant_pattern()`](https://wolski.github.io/prolfqua/reference/effective_contaminant_pattern.md),
[`effective_decoy_pattern()`](https://wolski.github.io/prolfqua/reference/effective_decoy_pattern.md),
[`get_uniprot_id_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
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
is_contaminant(c("zz|Cont00001|X", "sp|P2|X", "CON__ALBU"))
#> [1]  TRUE FALSE  TRUE
```
