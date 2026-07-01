# table facade to easily switch implementations

table facade to easily switch implementations

## Usage

``` r
table_facade(df, caption, digits = getOption("digits"), kable = TRUE)
```

## Arguments

- df:

  data.frame to display

- caption:

  table caption

- digits:

  number of digits (default from options)

- kable:

  if TRUE use knitr::kable

## Value

The requested plot, table, or transformed object.

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`effective_contaminant_pattern()`](https://wolski.github.io/prolfqua/reference/effective_contaminant_pattern.md),
[`effective_decoy_pattern()`](https://wolski.github.io/prolfqua/reference/effective_decoy_pattern.md),
[`get_uniprot_id_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`is_contaminant()`](https://wolski.github.io/prolfqua/reference/is_contaminant.md),
[`is_decoy()`](https://wolski.github.io/prolfqua/reference/is_decoy.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`remove_na_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md)

## Examples

``` r
table_facade(head(iris, 3), caption = "Iris preview")
#> 
#> 
#> Table: Iris preview
#> 
#> | Sepal.Length| Sepal.Width| Petal.Length| Petal.Width|Species |
#> |------------:|-----------:|------------:|-----------:|:-------|
#> |          5.1|         3.5|          1.4|         0.2|setosa  |
#> |          4.9|         3.0|          1.4|         0.2|setosa  |
#> |          4.7|         3.2|          1.3|         0.2|setosa  |
```
