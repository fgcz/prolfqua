# Detect decoy / reverse-database identifiers

Flags identifiers that look like decoy (reversed) database entries.
Built-in anchored default prefixes are always considered, unioned with
an optional configured `pattern`. An empty string or `NULL` pattern uses
the defaults only.

## Usage

``` r
is_decoy(ids, pattern = NULL)
```

## Arguments

- ids:

  character vector of (prefixed) identifiers

- pattern:

  optional additional decoy regex, unioned with the defaults

## Value

logical vector, `TRUE` where an id looks like a decoy

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`effective_contaminant_pattern()`](https://wolski.github.io/prolfqua/reference/effective_contaminant_pattern.md),
[`effective_decoy_pattern()`](https://wolski.github.io/prolfqua/reference/effective_decoy_pattern.md),
[`get_uniprot_id_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`is_contaminant()`](https://wolski.github.io/prolfqua/reference/is_contaminant.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`remove_na_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
is_decoy(c("REV_sp|P1|X", "sp|P2|X", "decoy_3", "normalProtein"))
#> [1]  TRUE FALSE  TRUE FALSE
is_decoy(c("shuffled_1", "P2"), pattern = "^shuffled_")
#> [1]  TRUE FALSE
```
