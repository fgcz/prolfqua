# generates peptide level plots for all Proteins

generates peptide level plots for all Proteins

## Usage

``` r
plot_hierarchies_boxplot_df(
  pdata,
  config,
  hierarchy = config$hierarchy_keys_depth(),
  facet_grid_on = NULL
)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

- hierarchy:

  e.g. protein_Id default hierarchy_keys_depth

- facet_grid_on:

  default NULL

## See also

Other plotting:
[`ContrastsPlotter`](https://wolski.github.io/prolfqua/reference/ContrastsPlotter.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`missigness_histogram()`](https://wolski.github.io/prolfqua/reference/missigness_histogram.md),
[`missingness_per_condition()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition_cumsum.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_heatmap()`](https://wolski.github.io/prolfqua/reference/plot_heatmap.md),
[`plot_heatmap_cor()`](https://wolski.github.io/prolfqua/reference/plot_heatmap_cor.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`plot_intensity_distribution_violin()`](https://wolski.github.io/prolfqua/reference/plot_intensity_distribution_violin.md),
[`plot_na_heatmap()`](https://wolski.github.io/prolfqua/reference/plot_NA_heatmap.md),
[`plot_pca()`](https://wolski.github.io/prolfqua/reference/plot_pca.md),
[`plot_raster()`](https://wolski.github.io/prolfqua/reference/plot_raster.md),
[`plot_sample_correlation()`](https://wolski.github.io/prolfqua/reference/plot_sample_correlation.md),
[`upset_interaction_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_interaction_missing_stats.md),
[`upset_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_missing_stats.md)

## Examples

``` r
 istar <- sim_lfq_data_peptide_config(with_missing = FALSE)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
 res <- plot_hierarchies_boxplot_df(istar$data,istar$config)
 istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
 config <- istar$config
 analysis <- istar$data
 analysis <- analysis |>
   dplyr::filter(protein_Id %in% sample(protein_Id, 2))

 res <- plot_hierarchies_boxplot_df(analysis,config)
 res$boxplot[[1]]
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 8 rows containing missing values or values outside the scale range
#> (`position_quasirandom()`).

 res <- plot_hierarchies_boxplot_df(analysis,config,config$hierarchy_keys()[1])
 res$boxplot[[1]]
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 8 rows containing missing values or values outside the scale range
#> (`position_quasirandom()`).

 res <- plot_hierarchies_boxplot_df(analysis,config,
                                    config$hierarchy_keys()[1],
                                    facet_grid_on = config$hierarchy_keys()[2])
 res$boxplot[[1]]
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 8 rows containing non-finite outside the scale range
#> (`stat_summary()`).
#> Warning: Removed 8 rows containing missing values or values outside the scale range
#> (`position_quasirandom()`).

 res$boxplot[[2]]


 iostar <- sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
 iostar$data <- iostar$data |>
   dplyr::filter(protein_Id %in% sample(protein_Id, 4))
 unique(iostar$data$protein_Id)
#> [1] "BEJI92~5282" "CGzoYe~2147" "SGIVBl~5782"

 res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config)
 res$boxplot[[1]]
#> Warning: Removed 1 row containing non-finite outside the scale range (`stat_boxplot()`).
#> Warning: Removed 1 row containing non-finite outside the scale range (`stat_summary()`).
#> Warning: Removed 1 row containing non-finite outside the scale range (`stat_summary()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`position_quasirandom()`).

 res <- plot_hierarchies_boxplot_df(iostar$data,iostar$config,
                                    iostar$config$hierarchy_keys()[1])
```
