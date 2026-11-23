# Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.

Compute mean, sd, and CV for all Peptides, or proteins, for all
interactions and all samples.

summarize stats output (compute quantiles)

## Usage

``` r
summarize_stats(lfqdata, factor_key = lfqdata$relevant_factor_keys())

summarize_stats_quantiles(
  stats_res,
  factor_keys_depth,
  stats = c("sd", "CV"),
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- factor_key:

  character vector - factor columns to group by (default:
  relevant_factor_keys)

- stats_res:

  result of running `summarize_stats`

- factor_keys_depth:

  character vector - factor columns at current depth

- stats:

  summarize either sd or CV

- probs:

  for which quantiles 10, 20 etc.

## Value

`summarize_stats()`: a tibble with one row per hierarchy element and
level of `factor_key`, with the columns `nrReplicates`, `nrMeasured`,
`nrNAs`, `sd`, `var`, `meanAbundance`, `medianAbundance`, `CV`
(untransformed data only) and `interaction`.

`summarize_stats_quantiles()`: a list with `long`, a tibble with one row
per factor level and quantile (columns `probs` and `quantiles`), and
`wide`, the same quantiles with one column per factor level.

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Stats()$stats()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Stats()$stats_quantiles()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
