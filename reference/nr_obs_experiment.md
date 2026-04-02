# Aggregates e.g. protein abundances from peptide abundances

Aggregates e.g. protein abundances from peptide abundances

## Usage

``` r
nr_obs_experiment(
  data,
  config,
  from_children = TRUE,
  name_nr_child = "nr_child_exp"
)
```

## Arguments

- data:

  tidy data

- config:

  prolfqua config

- from_children:

  TRUE compute from nr_children column, FALSE - do not use nr_children
  column

- name_nr_child:

  how to name column

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done


xd <- nr_obs_experiment(dd$data, dd$config, from_children = FALSE)
xd <- nr_obs_experiment(dd$data, dd$config, from_children = TRUE)
stopifnot(min(xd$nr_child_exp) == 1)
dp <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
nr_obs_experiment(dp$data, dp$config)
#> # A tibble: 10 × 2
#>    protein_Id  nr_child_exp
#>    <chr>              <dbl>
#>  1 0EfVhX~0087            3
#>  2 7cbcrd~5725            1
#>  3 9VUkAq~4703            1
#>  4 BEJI92~5282            2
#>  5 CGzoYe~2147            1
#>  6 DoWup2~5896            1
#>  7 Fl4JiV~8625            4
#>  8 HvIpHG~9079            2
#>  9 JcKVfU~9653            7
#> 10 SGIVBl~5782            6
nr_obs_experiment(dp$data, dp$config, from_children = FALSE)
#> # A tibble: 10 × 2
#>    protein_Id  nr_child_exp
#>    <chr>              <int>
#>  1 0EfVhX~0087            1
#>  2 7cbcrd~5725            1
#>  3 9VUkAq~4703            1
#>  4 BEJI92~5282            1
#>  5 CGzoYe~2147            1
#>  6 DoWup2~5896            1
#>  7 Fl4JiV~8625            1
#>  8 HvIpHG~9079            1
#>  9 JcKVfU~9653            1
#> 10 SGIVBl~5782            1


dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
dd$config$hierarchy_depth <- 2
xpep <- nr_obs_experiment(dd$data,dd$config)
stopifnot(all(xpep$nr_child_exp == 1))
xpep <- nr_obs_experiment(dd$data,dd$config, from_children = FALSE)
stopifnot(all(xpep$nr_child_exp == 1))
```
