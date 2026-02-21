# prepare benchmark data

prepare benchmark data

## Usage

``` r
ionstar_bench_preprocess(data, idcol = "protein_Id")
```

## Arguments

- data:

  analysis results

- idcol:

  default "protein_Id"

## See also

Other benchmarking:
[`Benchmark`](https://wolski.github.io/prolfqua/reference/Benchmark.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`make_benchmark()`](https://wolski.github.io/prolfqua/reference/make_benchmark.md),
[`ms_bench_add_scores()`](https://wolski.github.io/prolfqua/reference/ms_bench_add_scores.md),
[`ms_bench_auc()`](https://wolski.github.io/prolfqua/reference/ms_bench_auc.md)

## Examples

``` r
dd <- data.frame(
  protein_Id = c("P1_HUMAN", "P2_ECOLI", "P3_HUMAN", "P4_OTHER"),
  estimate = c(0.5, -1.2, 0.3, 0.1),
  stringsAsFactors = FALSE)
res <- ionstar_bench_preprocess(dd)
stopifnot(is.list(res))
stopifnot(all(c("data", "table") %in% names(res)))
stopifnot(nrow(res$data) == 3)  # OTHER is filtered out
stopifnot("TP" %in% colnames(res$data))
```
