# Compute theoretical sample sizes from factor level standard deviations

Compute theoretical sample sizes from factor level standard deviations

## Usage

``` r
lfq_power_t_test_proteins(
  stats_res,
  delta = c(0.59, 1, 2),
  power = 0.8,
  sig.level = 0.05,
  min.n = 1.5
)
```

## Arguments

- stats_res:

  data.frame `summarize_stats` output

- delta:

  effect size you are interested in

- power:

  of test

- sig.level:

  P-Value

- min.n:

  smallest n to determine

## Value

a tibble with one row per complete row of `stats_res` and value of
`delta`, with the added columns `delta`, `N_exact` (samples per group
needed) and `N` (rounded up)

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$get_Stats()$power_t_test()`.
See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
