# Format a table for reports

Formats a table for reports with
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html), or returns
it unchanged, so report code does not depend on the table package.

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

  if TRUE use knitr::kable (requires knitr)

## Value

a `knitr_kable` table, or `df` unchanged when `kable = FALSE`

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_uniprot_id_from_fasta_header()`](https://fgcz.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`is_contaminant()`](https://fgcz.github.io/prolfqua/reference/is_contaminant.md),
[`is_decoy()`](https://fgcz.github.io/prolfqua/reference/is_decoy.md)

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
