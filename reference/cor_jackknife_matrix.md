# Compute correlation matrix with jackknife resampling

Compute correlation matrix with jackknife resampling

## Usage

``` r
cor_jackknife_matrix(
  dataX,
  distmethod = function(x) {
     cor(x, use = "pairwise.complete.obs", method =
    "pearson")
 }
)
```

## Arguments

- dataX:

  data.frame with e.g. transition intensities per peptide or pepitde
  intensities per protein

## See also

Other transitioncorrlation:
[`cor_order()`](https://wolski.github.io/prolfqua/reference/cor_order.md),
[`jackknife()`](https://wolski.github.io/prolfqua/reference/jackknife.md),
[`jackknife_matrix()`](https://wolski.github.io/prolfqua/reference/jackknife_matrix.md)

## Examples

``` r
dd <- prolfqua_data("data_correlatedPeptideList")

class(dd[[1]])
#> [1] "data.frame"
dd[[1]][1,2] <- NA
cor_jackknife_matrix(dd[[1]])
#>                                             sp|O95477|ABCA1~HUMAN_EGAFVELFHEIDDR_3
#> sp|O95477|ABCA1~HUMAN_EGAFVELFHEIDDR_3                                   1.0000000
#> sp|O95477|ABCA1~HUMAN_FVSPLSWDLVGR_2                                     0.8674461
#> sp|O95477|ABCA1~HUMAN_LVEDIGHELTYVLPYEAAK_3                              0.9342006
#>                                             sp|O95477|ABCA1~HUMAN_FVSPLSWDLVGR_2
#> sp|O95477|ABCA1~HUMAN_EGAFVELFHEIDDR_3                                 0.8674461
#> sp|O95477|ABCA1~HUMAN_FVSPLSWDLVGR_2                                   1.0000000
#> sp|O95477|ABCA1~HUMAN_LVEDIGHELTYVLPYEAAK_3                            0.7940341
#>                                             sp|O95477|ABCA1~HUMAN_LVEDIGHELTYVLPYEAAK_3
#> sp|O95477|ABCA1~HUMAN_EGAFVELFHEIDDR_3                                        0.9342006
#> sp|O95477|ABCA1~HUMAN_FVSPLSWDLVGR_2                                          0.7940341
#> sp|O95477|ABCA1~HUMAN_LVEDIGHELTYVLPYEAAK_3                                   1.0000000
```
