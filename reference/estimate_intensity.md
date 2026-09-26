# Aggregates e.g. protein abundances from peptide abundances

Aggregates e.g. protein abundances from peptide abundances

## Usage

``` r
estimate_intensity(lfqdata, method = c("medpolish", "rlm"))
```

## Arguments

- lfqdata:

  LFQData object

- method:

  "medpolish" (Tukey's median polish, column \`medpolish\`) or "rlm"
  (robust regression with \`MASS::rlm\`, column \`lmrob\`)

## Value

returns list with data (data.frame) and config (AnalysisConfiguration)

## See also

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_top_n()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_top_n.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(dd$data, dd$config)
lfq <- lfq$get_Transformer()$log2()$lfq
#> Column added : log2_abundance
bbMed <- estimate_intensity(lfq, method = "medpolish")
#> starting aggregation
bbRob <- estimate_intensity(lfq, method = "rlm")
#> starting aggregation
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
nrow(bbMed$data)
#> [1] 116
nrow(bbRob$data)
#> [1] 116
xt <- dplyr::inner_join(bbMed$data, bbRob$data)
#> Joining with `by = join_by(protein_Id, sampleName, group_, sample,
#> isotopeLabel, nr_children_protein_Id)`
plot(xt$medpolish, xt$lmrob, log = "xy", pch = "*")
abline(0, 1, col = 2)

```
