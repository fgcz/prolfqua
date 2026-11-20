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
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, group, mean, N, Replicate.
#> setup_analysis: added 22 rows with NA response for feature and sample combinations absent from the input.
res2$config$factor_depth <- 2
lfq2 <- LFQData$new(res2$data, res2$config)
xx <- summarize_stats_factors(lfq2)
stopifnot(nrow(xx) == 80)
stopifnot(length(unique(xx$interaction)) == (2 + 2 + 2 * 2))
```
