# Aggregates top N intensities

run
[rank_peptide_by_intensity](https://fgcz.github.io/prolfqua/reference/rank_peptide_by_intensity.md)
first

## Usage

``` r
aggregate_intensity_top_n(ranked_data, lfqdata, .func, N = 3)
```

## Arguments

- ranked_data:

  data.frame with ranked peptides

- lfqdata:

  LFQData object

- .func:

  function to use for aggregation

- N:

  default 3 top intensities.

## Value

list with data and new reduced configuration (config)

## See also

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`estimate_intensity()`](https://fgcz.github.io/prolfqua/reference/estimate_intensity.md),
[`plot_estimate()`](https://fgcz.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_line()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_line.md)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(dd$data, dd$config)
ranked <- rank_peptide_by_intensity(lfq$data_long(), lfq$response(), lfq$hierarchy_keys())
#> Columns added : srm_meanInt srm_meanIntRank

mean_f <- function(x, name = FALSE) {
  if (name) return("mean")
  mean(x, na.rm = TRUE)
}

resTOPN <- aggregate_intensity_top_n(ranked, lfq, .func = mean_f, N = 3)
stopifnot(names(resTOPN) %in% c("data", "config"))
lfq_agg <- LFQData$new(resTOPN$data, resTOPN$config)
tmpRob <- plot_estimate(lfq, lfq_agg, show.legend = TRUE)
stopifnot("ggplot" %in% class(tmpRob$plots[[4]]))
```
