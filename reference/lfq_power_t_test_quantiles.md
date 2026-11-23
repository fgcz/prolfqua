# Estimate sample sizes from standard deviation quantiles

Estimate sample sizes from standard deviation quantiles

## Usage

``` r
lfq_power_t_test_quantiles(
  quantile_sd,
  delta = c(0.59, 1, 2),
  power = 0.8,
  sig.level = 0.05,
  min.n = 1.5
)
```

## Arguments

- quantile_sd:

  output of `summarize_stats_quantiles`

- delta:

  effect size you are interested in

- power:

  of test

- sig.level:

  P-Value

- min.n:

  smallest n to determine

## Value

a tibble with the rows of `quantile_sd` repeated for each `delta`, with
the added columns `sdtrimmed` (the sd floored at the smallest sd
resolvable with `min.n` samples), `N_exact`, `N` (rounded up) and
`delta`

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Stats()$power_t_test_quantiles()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
