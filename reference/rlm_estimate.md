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

  e.g. sampleName

## See also

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_topN()`](https://wolski.github.io/prolfqua/reference/aggregate_intensity_topN.md),
[`estimate_intensity()`](https://wolski.github.io/prolfqua/reference/estimate_intensity.md),
[`intensity_summary_by_hkeys()`](https://wolski.github.io/prolfqua/reference/intensity_summary_by_hkeys.md),
[`medpolish_estimate()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate.md),
[`medpolish_estimate_df()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_df.md),
[`medpolish_estimate_dfconfig()`](https://wolski.github.io/prolfqua/reference/medpolish_estimate_dfconfig.md),
[`medpolish_protein_estimates()`](https://wolski.github.io/prolfqua/reference/medpolish_protein_estimates.md),
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
  sampleName = rep(letters[1:4], 5))
rlm_estimate(xx2, "log2Area", "peptide_Id", "sampleName")
#>   sampleName mean.log2Area    weights     lmrob
#> 1          a     -3.562204 0.01693638 -3.198599
#> 2          b      2.928612 0.02229354  2.928612
#> 3          c     -1.869415 0.02596861 -1.869415
#> 4          d     -2.967742 0.01037056 -4.903623
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
  res[[i]] <- list(data = xx3, summary = rlm_estimate(xx3, "log2Area", "peptide_Id", "sampleName"))
}
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
#> Warning: 'rlm' failed to converge in 20 steps
rlm_estimate(xx2[xx2$peptide_Id == "A", ], "log2Area", "peptide_Id", "sampleName")
#>    sampleName mean.log2Area       lmrob weights
#> 1           a   -21.3168544 -21.3168544       1
#> 6           b    -8.0899700  -8.0899700       1
#> 11          c    -0.9207308  -0.9207308       1
#> 16          d    11.0102264  11.0102264       1
rlm_estimate(xx2[xx2$sampleName == "a", ], "log2Area", "peptide_Id", "sampleName")
#> # A tibble: 1 × 4
#>   sampleName lmrob mean.log2Area weights
#>   <chr>      <dbl>         <dbl>   <dbl>
#> 1 a          -3.56         -3.56       1


bb <- prolfqua_data("data_ionstar")$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
bb$config <- old2new(bb$config)
stopifnot(nrow(bb$data) == 25780)
conf <- bb$config
data <- bb$data
conf$hierarchyDepth <- 1
xnested <- data |>
  dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
  tidyr::nest()

feature <- base::setdiff(conf$hierarchy_keys(), conf$hierarchy_keys_depth())
x <- xnested$data[[1]]
bb <- rlm_estimate(x,
  response = conf$get_response(),
  feature = feature,
  samples = conf$sampleName
)
```
