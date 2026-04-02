# Aggregates e.g. protein abundances from peptide abundances

Aggregates e.g. protein abundances from peptide abundances

## Usage

``` r
estimate_intensity(data, config, .func)
```

## Arguments

- func:

  \- a function working on a matrix of intensities for each protein.

## Value

returns list with data (data.frame) and config (AnalysisConfiguration)

## See also

[`medpolish_estimate_dfconfig`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md)
[`rlm_estimate_dfconfig`](https://wolski.github.io/prolfqua/reference/rlm_estimate_dfconfig.md)

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_topN()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_topN.md),
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
data <- dd$data
data <- prolfqua::transform_work_intensity(data, config, log2)
#> Column added : log2_abundance
bbMed <- estimate_intensity(data, config, .func = medpolish_estimate_dfconfig)
#> starting aggregation
bbRob <- estimate_intensity(data, config, .func = rlm_estimate_dfconfig)
#> starting aggregation
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
nrow(bbMed$data)
#> [1] 116
nrow(bbRob$data)
#> [1] 116
length(bbMed$data$medpolish)
#> [1] 116
length(bbRob$data$lmrob)
#> [1] 116
xt <- dplyr::inner_join(bbMed$data, bbRob$data)
#> Joining with `by = join_by(protein_Id, sampleName, group_, sample,
#> isotopeLabel, nr_children_protein_Id)`
plot(xt$medpolish, xt$lmrob, log = "xy", pch = "*")
abline(0, 1, col = 2)

```
