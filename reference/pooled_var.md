# pooled variance

following the documentation here:
https://online.stat.psu.edu/stat500/lesson/7/7.3/7.3.1/7.3.1.1

## Usage

``` r
poolvar(res1, config)

compute_pooled(x)
```

## Arguments

- res1:

  per-group statistics from
  [`summarize_stats()`](https://fgcz.github.io/prolfqua/reference/summarize_stats.md)

- config:

  AnalysisConfiguration of the data

- x:

  data.frame

## Value

data.frame

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Stats()$power_fdr_pilot()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Stats()$power_fdr_pilot()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
