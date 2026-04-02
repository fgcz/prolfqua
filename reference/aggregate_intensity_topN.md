# Aggregates top N intensities

run
[rank_peptide_by_intensity](https://wolski.github.io/prolfqua/reference/rank_peptide_by_intensity.md)
first

## Usage

``` r
aggregate_intensity_topN(pdata, config, .func, N = 3)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

- N:

  default 3 top intensities.

- func:

  function to use for aggregation

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
config <- dd$config
res <- dd$data
ranked <- rank_peptide_by_intensity(res, config)
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank

mean_f <- function(x, name = FALSE) {
  if (name) {
    return("mean")
  }
  mean(x, na.rm = TRUE)
}
sum_f <- function(x, name = FALSE) {
  if (name) {
    return("sum")
  }
  sum(x, na.rm = TRUE)
}

resTOPN <- aggregate_intensity_topN(
  ranked,
  config,
  .func = mean_f,
  N = 3
)

print(dim(resTOPN$data))
#> [1] 116   8
# stopifnot(dim(resTOPN$data) == c(3260, 8))
stopifnot(names(resTOPN) %in% c("data", "config"))
config$get_response()
#> [1] "abundance"
tmpRob <- plot_estimate(ranked,
  config,
  resTOPN$data,
  resTOPN$config,
  show.legend = TRUE
)
stopifnot("ggplot" %in% class(tmpRob$plots[[4]]))
```
