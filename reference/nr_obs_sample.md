# Aggregates e.g. protein abundances from peptide abundances

Aggregates e.g. protein abundances from peptide abundances

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

  character — intensity column name

- hierarchy_keys_depth:

  character vector — hierarchy columns at current depth

- file_name:

  character — file name column

- nr_children_col:

  character — nr_children column name

- new_child:

  character — output column name

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
dd$data <- na.omit(dd$data)

xd <- nr_obs_sample(na.omit(dd$data), dd$config$get_response(),
  dd$config$hierarchy_keys_depth(), dd$config$file_name, dd$config$nr_children)
xd$nr_children |> table()
#> 
#>  1  2  3  4  5  6  7 
#> 52 24 11  7  5  5 12 
# xd |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_children)

dp <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
xp <- nr_obs_sample(dp$data, dp$config$get_response(),
  dp$config$hierarchy_keys_depth(), dp$config$file_name, dp$config$nr_children)
# xp
# xp |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_peptides)
xp$nr_peptides |> table()
#> 
#>  1  2  3  4  6  7 
#> 37 22 11 11 11 12 
```
