# Aggregates e.g. protein abundances from peptide abundances

Aggregates e.g. protein abundances from peptide abundances

## Usage

``` r
nr_obs_sample(data, config, new_child = config$nr_children)
```

## Arguments

- data:

  data.frame with peptide-level data

- config:

  AnalysisConfiguration

- new_child:

  column name for the nr_children count

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
dd$data <- na.omit(dd$data)

xd <- nr_obs_sample(dd$data, dd$config)
xd$nr_children |> table()
#> 
#>  1  2  3  4  5  6  7 
#> 52 24 11  7  5  5 12 
# xd |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_children)

dp <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
xp <- nr_obs_sample(dp$data, dp$config)
# xp
# xp |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_peptides)
xp$nr_peptides |> table()
#> 
#>  1  2  3  4  6  7 
#> 37 22 11 11 11 12 
```
