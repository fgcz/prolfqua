# Keep only the proteins with enough peptides

Keeps the rows of the proteins with at least
\`config\$min_peptides_protein\` children at the next hierarchy level,
counted by \[nr_B_in_A()\], and adds that count. The package dataset
\`data_ionstar\` calls it from its \`filtered()\` method.

## Usage

``` r
filter_proteins_by_peptide_count(pdata, config)
```

## Arguments

- pdata:

  data.frame in the long format of \`config\`

- config:

  AnalysisConfiguration

## Value

list with the filtered \`data\` and the \`name\` of the count column,
\`NULL\` when the hierarchy has no child level

## See also

Other preprocessing:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`apply_to_response_matrix()`](https://fgcz.github.io/prolfqua/reference/apply_to_response_matrix.md),
[`get_robscales()`](https://fgcz.github.io/prolfqua/reference/get_robscales.md),
[`nr_B_in_A()`](https://fgcz.github.io/prolfqua/reference/nr_B_in_A.md),
[`robust_scale()`](https://fgcz.github.io/prolfqua/reference/robust_scale.md),
[`scale_with_subset()`](https://fgcz.github.io/prolfqua/reference/scale_with_subset.md)

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
filtered <- prolfqua:::filter_proteins_by_peptide_count(istar$data, istar$config)
stopifnot(all(filtered$data[[filtered$name]] >= istar$config$min_peptides_protein))
```
