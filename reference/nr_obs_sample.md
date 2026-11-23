# Sum of observed children per hierarchy unit and sample

For each hierarchy unit and sample, sums the `nr_children_col` column
over the rows with a non-missing response, e.g. the number of peptides
behind each protein intensity of a sample. Counts per sample, observed
values only.

## Usage

``` r
nr_obs_sample(
  data,
  response,
  hierarchy_keys_depth,
  file_name,
  nr_children_col,
  new_child = nr_children_col
)
```

## Arguments

- data:

  data.frame

- response:

  character - intensity column name

- hierarchy_keys_depth:

  character vector - hierarchy columns at current depth

- file_name:

  character - file name column

- nr_children_col:

  character - nr_children column name

- new_child:

  character - output column name

## Value

a tibble with one row per hierarchy element and sample:
`hierarchy_keys_depth`, `file_name` and `new_child`, the sum of
`nr_children_col` over rows with an observed response.

## See also

Other hierarchy counting:
[`HierarchyCountsSample`](https://fgcz.github.io/prolfqua/reference/HierarchyCountsSample.md),
[`nr_children_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_children_experiment.md),
[`nr_features_experiment()`](https://fgcz.github.io/prolfqua/reference/nr_features_experiment.md)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
dd$data <- na.omit(dd$data)

xd <- nr_obs_sample(na.omit(dd$data), dd$config$get_response(),
  dd$config$hierarchy_keys_depth(), dd$config$file_name, dd$config$nr_children)
xd$nr_children |> table()
#> 
#>  1  2  3  4  5  6  7 
#> 52 24 11  7  5  5 12 
# xd |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_children)

dp <- prolfqua::sim_lfq_data_protein_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 16 rows with NA response for feature and sample combinations absent from the input.
xp <- nr_obs_sample(dp$data, dp$config$get_response(),
  dp$config$hierarchy_keys_depth(), dp$config$file_name, dp$config$nr_children)
# xp
# xp |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nrPeptides)
xp$nrPeptides |> table()
#> 
#>  1  2  3  4  6  7 
#> 37 22 11 11 11 12 
```
