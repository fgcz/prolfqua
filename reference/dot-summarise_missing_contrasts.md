# summarise_missing_contrasts

summarise_missing_contrasts

## Usage

``` r
.summarise_missing_contrasts(
  data,
  hierarchy = c("protein_Id"),
  contrast = "contrast",
  what = "statistic"
)
```

## Examples

``` r
ttd <- ionstar_bench_preprocess(prolfqua_data('data_benchmarkExample'))
x <- .summarise_missing_contrasts(ttd$data)
x2 <- tibble::as_tibble(x$summary)
```
