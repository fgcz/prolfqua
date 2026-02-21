# adds FDR, TPR and FDP to data.

adds FDR, TPR and FDP to data.

## Usage

``` r
ms_bench_add_scores(
  data,
  TP_col = "TP",
  arrangeby = "diff",
  desc = TRUE,
  subject_Id = "protein_Id"
)
```

## Arguments

- data:

  a dataframe with TP_col indicating if the row is a true positive hit

- TP_col:

  column name of TP (TRUE, FALSE)

- arrangeby:

  \- by which column to sort.

- desc:

  descending or ascending.

## Value

data.frame with the following columns added FDP - false discovery
proportion (Q in Benjamini Hochberg table) FPR - false positive rate
TPR - true positive rate TP_hits - true positives

## See also

Other benchmarking:
[`Benchmark`](https://wolski.github.io/prolfqua/reference/Benchmark.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`ionstar_bench_preprocess()`](https://wolski.github.io/prolfqua/reference/ionstar_bench_preprocess.md),
[`make_benchmark()`](https://wolski.github.io/prolfqua/reference/make_benchmark.md),
[`ms_bench_auc()`](https://wolski.github.io/prolfqua/reference/ms_bench_auc.md)

## Examples

``` r
dd <- prolfqua_data('data_test_confusion_matrix_scores')
xd <- ms_bench_add_scores(dd, arrangeby = "estimate")
plot(xd$TPR,xd$PREC, type="l")

plot(1- xd$PREC, xd$FDP)

```
