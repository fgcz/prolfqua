# Extracts uniprot ID

Extracts uniprot ID

## Usage

``` r
get_UniprotID_from_fasta_header(df, idcolumn = "protein_Id")
```

## Arguments

- df:

  data.frame

- idcolumn:

  character, column name to extract uniprot IDs from

## Value

data.frame

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`remove_NA_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
bb <- prolfqua_data('data_ionstar')$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
tmp <- prolfqua::separate_hierarchy(bb$data,old2new( bb$config))
tmp$UniprotID <- NULL
tmp <- get_UniprotID_from_fasta_header(tmp, idcolumn = "top_protein")
stopifnot("UniprotID" %in%  colnames(tmp))
```
