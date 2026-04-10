# Estimate e.g. protein abundance from peptides using MASS:rlm

Estimate e.g. protein abundance from peptides using MASS:rlm

## Usage

``` r
rlm_estimate(pdata, response, feature, samples, maxIt = 20)
```

## Arguments

- pdata:

  data

- response:

  intensities

- feature:

  e.g. peptideIDs.

- samples:

  e.g. sample_name

## See also

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_top_n()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_top_n.md),
[`estimate_intensity()`](https://wolski.github.io/prolfqua/reference/estimate_intensity.md),
[`medpolish_estimate()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`medpolish_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md),
[`plot_estimate()`](https://wolski.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_add_quantline()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_add_quantline.md),
[`plot_hierarchies_line()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line.md),
[`plot_hierarchies_line_df()`](https://wolski.github.io/prolfqua/reference/plot_hierarchies_line_df.md),
[`rlm_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/rlm_estimate_dfconfig.md)

## Examples

``` r
xx <- data.frame(response = rnorm(20, 0, 10), feature = rep(LETTERS[1:5], 4),
  samples = rep(letters[1:4], 5))

bb <- rlm_estimate(xx, "response", "feature", "samples", maxIt = 20)

xx2 <- data.frame(log2Area = rnorm(20, 0, 10), peptide_Id = rep(LETTERS[1:5], 4),
  sample_name = rep(letters[1:4], 5))
rlm_estimate(xx2, "log2Area", "peptide_Id", "sample_name")
#>   sample_name mean.log2Area    weights     lmrob
#> 1           a    -0.5323711 0.01276912 -1.674944
#> 2           b    -3.5622039 0.01458104 -3.383001
#> 3           c     2.9286124 0.02616560  2.928612
#> 4           d    -1.8694154 0.02666118 -1.869415
rlm_estimate(prolfqua_data("data_checksummarizationrobust87"),
  "log2Area", "peptide_Id", "sampleName")
#>   sampleName mean.log2Area weights      lmrob
#> 1          a    -7.1898238       1 -7.1898238
#> 2          b    -6.3420890       1 -6.3420890
#> 3          c            NA      NA         NA
#> 4          d     0.4191734       1  0.4191734
rlm_estimate(prolfqua_data("data_checksummarizerobust69"),
  "log2Area", "peptide_Id", "sampleName")
#>   sampleName mean.log2Area    weights     lmrob
#> 1          a     -6.307444 0.05027647 -6.307444
#> 2          b            NA         NA        NA
#> 3          c      2.489945 0.00329409  2.489945
#> 4          d     -9.076788 0.01213768 -9.076788
res <- vector(100, mode = "list")
for (i in seq_len(100)) {
  xx3 <- xx2
  xx3$log2Area[sample(1:20, sample(1:15, 1))] <- NA
  res[[i]] <- list(data = xx3, summary = rlm_estimate(xx3, "log2Area", "peptide_Id", "sample_name"))
}
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
rlm_estimate(xx2[xx2$peptide_Id == "A", ], "log2Area", "peptide_Id", "sample_name")
#>    sample_name mean.log2Area       lmrob weights
#> 1            a    4.06124230  4.06124230       1
#> 6            b   -5.76691823 -5.76691823       1
#> 11           c   12.40584606 12.40584606       1
#> 16           d   -0.01408112 -0.01408112       1
rlm_estimate(xx2[xx2$sample_name == "a", ], "log2Area", "peptide_Id", "sample_name")
#> # A tibble: 1 × 4
#>   sample_name  lmrob mean.log2Area weights
#>   <chr>        <dbl>         <dbl>   <dbl>
#> 1 a           -0.532        -0.532       1


bb <- sim_lfq_data_peptide_config(Nprot = 20)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
conf <- bb$config
data <- bb$data
conf$hierarchy_depth <- 1
xnested <- data |>
  dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
  tidyr::nest()

feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
x <- xnested$data[[1]]
bb <- rlm_estimate(x,
  response = conf$get_response(),
  feature = feature,
  samples = conf$sample_name
)
```
