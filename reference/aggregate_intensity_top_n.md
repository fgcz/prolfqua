# Aggregates top N intensities

run
[rank_peptide_by_intensity](https://wolski.github.io/prolfqua/reference/rank_peptide_by_intensity.md)
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
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`estimate_intensity()`](https://wolski.github.io/prolfqua/reference/estimate_intensity.md),
[`medpolish_estimate()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`medpolish_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`rlm_estimate()`](https://wolski.github.io/prolfqua/reference/rlm_estimate.md),
[`rlm_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/rlm_estimate_dfconfig.md)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(dd$data, dd$config)
ranked <- rank_peptide_by_intensity(lfq$get_data(), lfq$response(), lfq$hierarchy_keys())
#> Joining with `by = join_by(protein_Id, peptide_Id)`
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
