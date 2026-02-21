# make Benchmark

make Benchmark

## Usage

``` r
make_benchmark(
  prpr,
  contrast = "contrast",
  toscale = c("p.value"),
  avgInt = "avgInt",
  fcestimate = "diff",
  benchmark = list(list(score = "diff", desc = TRUE), list(score = "statistic", desc =
    TRUE), list(score = "scaled.p.value", desc = TRUE)),
  FDRvsFDP = list(list(score = "FDR", desc = FALSE)),
  model_description = "protein level measurments, linear model",
  model_name = "prot_med_lm",
  hierarchy = c("protein_Id"),
  summarizeNA = "statistic"
)
```

## Arguments

- prpr:

  prepared data, e.g. function
  [`ionstar_bench_preprocess`](https://wolski.github.io/prolfqua/reference/ionstar_bench_preprocess.md)

- contrast:

  column with names of the contrast

- toscale:

  which scores to scale using fcestimate, typically p.value

- avgInt:

  column with average intensity

- fcestimate:

  column with FC estimate

- benchmark:

  which scores to benchmark e.g. diff, statistics

- FDRvsFDP:

  which score to plot against the false discovery proportion

- model_description:

  string describing the model

- model_name:

  name of the model, compatible with
  [`make.names`](https://rdrr.io/r/base/make.names.html) output.

- hierarchy:

  name of column with protein ID.

- summarizeNA:

  summarizeNA

## Value

Benchmark

## See also

Other benchmarking:
[`Benchmark`](https://wolski.github.io/prolfqua/reference/Benchmark.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`ionstar_bench_preprocess()`](https://wolski.github.io/prolfqua/reference/ionstar_bench_preprocess.md),
[`ms_bench_add_scores()`](https://wolski.github.io/prolfqua/reference/ms_bench_add_scores.md),
[`ms_bench_auc()`](https://wolski.github.io/prolfqua/reference/ms_bench_auc.md)
