# Keep only the parents with enough children

Keeps the rows of the \`parent\` entries (e.g. proteins) with at least
\`threshold\` distinct \`child\` entries (e.g. peptides), counted by
\[nr_B_in_A()\], and adds that count. By default it counts the second
hierarchy level within the first, independent of \`hierarchy_depth\`.
The package dataset \`data_ionstar\` calls it as \`(pdata, config)\`
from its \`filtered()\` method.

## Usage

``` r
filter_proteins_by_peptide_count(
  pdata,
  config,
  parent = config$hierarchy_keys()[1],
  child = config$hierarchy_keys()[2],
  threshold = config$min_peptides_protein
)
```

## Arguments

- pdata:

  data.frame in the long format of \`config\`

- config:

  AnalysisConfiguration

- parent:

  parent key column, default the first hierarchy key

- child:

  child key column, default the second hierarchy key

- threshold:

  minimum number of children, default \`config\$min_peptides_protein\`

## Value

list with the filtered \`data\` and the \`name\` of the count column

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
filtered <- prolfqua:::filter_proteins_by_peptide_count(
  istar$data, istar$config, parent = "protein_Id", child = "peptide_Id", threshold = 2
)
stopifnot(all(filtered$data[[filtered$name]] >= 2))
```
