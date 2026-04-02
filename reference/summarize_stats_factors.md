# compute var sd etc for all factor levels

compute var sd etc for all factor levels

## Usage

``` r
summarize_stats_factors(pdata, config)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

## Examples

``` r
# example code
res2 <- prolfqua::sim_lfq_data_2Factor_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
xx <- summarize_stats_factors(res2$data, res2$config)
#> completing cases
stopifnot(nrow(xx) == 80)
stopifnot( length(unique(xx$interaction)) == (2 + 2 + 2 * 2))
```
