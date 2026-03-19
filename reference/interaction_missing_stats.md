# compute missingness statistics per hierarchy and factor level

compute missingness statistics per hierarchy and factor level

## Usage

``` r
interaction_missing_stats(
  pdata,
  config,
  factors = config$factor_keys_depth(),
  hierarchy = config$hierarchy_keys(),
  workIntensity = config$get_response()
)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

- factors:

  factor to include (default up to factor depth)

- hierarchy:

  hierarchy to include (default up to hierarchy depth)

- workIntensity:

  work intensity column

## Examples

``` r

istar <- sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- istar$config
analysis <- istar$data

xx <- complete_cases(analysis, config)
#> completing cases
x <- interaction_missing_stats(xx, config)$data |> dplyr::arrange(desc(nrNAs))
#> Warning: >>>> deprecated! <<<< 
#> 
#>           use summarize_stats_factors instead.
#> completing cases
nrow(x)
#> [1] 84
tmp <- interaction_missing_stats(xx, config,
 factors= character(),
  hierarchy = config$hierarchy_keys()[1])$data
#> Warning: >>>> deprecated! <<<< 
#> 
#>           use summarize_stats_factors instead.
#> completing cases
stopifnot(nrow(tmp) == 10)
tmp <- interaction_missing_stats(xx, config,
  hierarchy = config$hierarchy_keys()[1])$data
#> Warning: >>>> deprecated! <<<< 
#> 
#>           use summarize_stats_factors instead.
#> completing cases
stopifnot(nrow(tmp) == length(unique(xx$protein_Id))* length(unique(xx$group_)))
stopifnot(sum(is.na(tmp$nrMeasured))==0)

tmp <- interaction_missing_stats(xx, config, factors = NULL)
#> Warning: >>>> deprecated! <<<< 
#> 
#>           use summarize_stats_factors instead.
#> completing cases
```
