# computes auc and pauc using trapez rule

computes auc and pauc using trapez rule

## Usage

``` r
ms_bench_auc(FPR, TPR, fpr_threshold = 1)
```

## Arguments

- FPR:

  array of FPR

- TPR:

  array of corresponding TPR

- fpr_threshold:

  default = 1

## See also

Other benchmarking:
[`Benchmark`](https://wolski.github.io/prolfqua/reference/Benchmark.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`ionstar_bench_preprocess()`](https://wolski.github.io/prolfqua/reference/ionstar_bench_preprocess.md),
[`make_benchmark()`](https://wolski.github.io/prolfqua/reference/make_benchmark.md),
[`ms_bench_add_scores()`](https://wolski.github.io/prolfqua/reference/ms_bench_add_scores.md)

## Examples

``` r
FPR <- c(0, 0.1, 0.2, 0.5, 1.0)
TPR <- c(0, 0.4, 0.7, 0.9, 1.0)
auc <- ms_bench_auc(FPR, TPR)
stopifnot(is.numeric(auc), length(auc) == 1, auc > 0)
```
