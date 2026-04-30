# compute var sd etc for all factor levels

compute var sd etc for all factor levels

## Usage

``` r
summarize_stats_factors(lfqdata)
```

## Arguments

- lfqdata:

  LFQData object

## Value

The computed result.

## Examples

``` r
res2 <- prolfqua::sim_lfq_data_2factor_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
res2$config$factor_depth <- 2
lfq2 <- LFQData$new(res2$data, res2$config)
xx <- summarize_stats_factors(lfq2)
stopifnot(nrow(xx) == 80)
stopifnot(length(unique(xx$interaction)) == (2 + 2 + 2 * 2))
```
