# Compute correlation matrix

Compute correlation matrix

## Usage

``` r
cor_order(dataX)
```

## Arguments

- dataX:

  data.frame with transition intensities per peptide

## See also

Other transitioncorrlation:
[`cor_jackknife_matrix()`](https://wolski.github.io/prolfqua/reference/cor_jackknife_matrix.md),
[`jackknife()`](https://wolski.github.io/prolfqua/reference/jackknife.md),
[`jackknife_matrix()`](https://wolski.github.io/prolfqua/reference/jackknife_matrix.md)

## Examples

``` r
list <- prolfqua_data('data_correlatedPeptideList')
cor_order(data_correlatedPeptideList[[1]])
#>                                             sp|O95477|ABCA1~HUMAN_FVSPLSWDLVGR_2
#> sp|O95477|ABCA1~HUMAN_FVSPLSWDLVGR_2                                   1.0000000
#> sp|O95477|ABCA1~HUMAN_EGAFVELFHEIDDR_3                                 0.6352343
#> sp|O95477|ABCA1~HUMAN_LVEDIGHELTYVLPYEAAK_3                            0.6950875
#>                                             sp|O95477|ABCA1~HUMAN_EGAFVELFHEIDDR_3
#> sp|O95477|ABCA1~HUMAN_FVSPLSWDLVGR_2                                     0.6352343
#> sp|O95477|ABCA1~HUMAN_EGAFVELFHEIDDR_3                                   1.0000000
#> sp|O95477|ABCA1~HUMAN_LVEDIGHELTYVLPYEAAK_3                              0.8983625
#>                                             sp|O95477|ABCA1~HUMAN_LVEDIGHELTYVLPYEAAK_3
#> sp|O95477|ABCA1~HUMAN_FVSPLSWDLVGR_2                                          0.6950875
#> sp|O95477|ABCA1~HUMAN_EGAFVELFHEIDDR_3                                        0.8983625
#> sp|O95477|ABCA1~HUMAN_LVEDIGHELTYVLPYEAAK_3                                   1.0000000
```
