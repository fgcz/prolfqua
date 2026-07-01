# converts LFQData object to SummarizedExperiment

For compatibility with Bioconductor

## Usage

``` r
LFQDataToSummarizedExperiment(lfqdata)
```

## Arguments

- lfqdata:

  LFQData object

## Value

SummarizedExperiment (bioconductor)

## See also

Other LFQData:
[`AggregateLimpa`](https://wolski.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://wolski.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://wolski.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://wolski.github.io/prolfqua/reference/AggregateTopN.md),
[`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://wolski.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://wolski.github.io/prolfqua/reference/LFQDataSummariser.md)

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
data <- istar$data
lfqdata <- LFQData$new(data, istar$config)
lfqdata$data_wide()
#> $data
#> # A tibble: 28 × 15
#>    protein_Id  peptide_Id isotopeLabel  A_V1  A_V2  A_V3  A_V4  B_V1  B_V2  B_V3
#>    <chr>       <chr>      <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 0EfVhX~0087 ITLb4x1q   light         19.4  NA    18.2  18.9  28.9  27.3  29.9
#>  2 0EfVhX~0087 ahQLlQY7   light         26.9  25.2  26.1  25.2  26.1  26.2  24.5
#>  3 0EfVhX~0087 dJkdz7so   light         NA    16.1  14.7  NA    25.7  23.9  25.5
#>  4 7cbcrd~5725 D5dQ4nKk   light         30.4  30.4  29.6  27.8  NA    21.6  20.9
#>  5 9VUkAq~4703 eIC06D7g   light         17.1  17.3  19.0  18.8  16.4  NA    17.4
#>  6 BEJI92~5282 HBkZvdhT   light         18.1  17.9  NA    NA    19.1  19.2  18.7
#>  7 BEJI92~5282 qQ1GK8Un   light         21.7  24.0  24.7  24.8  28.6  25.8  26.2
#>  8 CGzoYe~2147 mjHSHhoe   light         24.1  24.2  25.0  23.8  NA    31.3  31.1
#>  9 DoWup2~5896 KVUnZ6oZ   light         24.2  24.4  23.3  23.4  18.4  18.5  18.4
#> 10 Fl4JiV~8625 GsUIOl6Q   light         20.2  20.1  20.3  NA    25.0  25.5  24.9
#> # ℹ 18 more rows
#> # ℹ 5 more variables: B_V4 <dbl>, Ctrl_V1 <dbl>, Ctrl_V2 <dbl>, Ctrl_V3 <dbl>,
#> #   Ctrl_V4 <dbl>
#> 
#> $annotation
#> # A tibble: 12 × 4
#>    sampleName sample  group_ isotopeLabel
#>    <chr>      <chr>   <chr>  <chr>       
#>  1 A_V1       A_V1    A      light       
#>  2 A_V2       A_V2    A      light       
#>  3 A_V3       A_V3    A      light       
#>  4 A_V4       A_V4    A      light       
#>  5 B_V1       B_V1    B      light       
#>  6 B_V2       B_V2    B      light       
#>  7 B_V3       B_V3    B      light       
#>  8 B_V4       B_V4    B      light       
#>  9 Ctrl_V1    Ctrl_V1 Ctrl   light       
#> 10 Ctrl_V2    Ctrl_V2 Ctrl   light       
#> 11 Ctrl_V3    Ctrl_V3 Ctrl   light       
#> 12 Ctrl_V4    Ctrl_V4 Ctrl   light       
#> 
#> $rowdata
#> # A tibble: 28 × 3
#>    protein_Id  peptide_Id isotopeLabel
#>    <chr>       <chr>      <chr>       
#>  1 0EfVhX~0087 ITLb4x1q   light       
#>  2 0EfVhX~0087 ahQLlQY7   light       
#>  3 0EfVhX~0087 dJkdz7so   light       
#>  4 7cbcrd~5725 D5dQ4nKk   light       
#>  5 9VUkAq~4703 eIC06D7g   light       
#>  6 BEJI92~5282 HBkZvdhT   light       
#>  7 BEJI92~5282 qQ1GK8Un   light       
#>  8 CGzoYe~2147 mjHSHhoe   light       
#>  9 DoWup2~5896 KVUnZ6oZ   light       
#> 10 Fl4JiV~8625 GsUIOl6Q   light       
#> # ℹ 18 more rows
#> 
#> $config
#> <AnalysisConfiguration>
#>   Public:
#>     annotation_vars: function () 
#>     bin_resp: 
#>     clone: function (deep = FALSE) 
#>     factor_depth: 1
#>     factor_keys: function () 
#>     factor_keys_depth: function () 
#>     factors: list
#>     file_name: sample
#>     get_response: function () 
#>     hierarchy: list
#>     hierarchy_depth: 1
#>     hierarchy_keys: function (rev = FALSE) 
#>     hierarchy_keys_depth: function (names = TRUE) 
#>     id_required: function () 
#>     id_vars: function () 
#>     ident_q_value: qValue
#>     ident_score: 
#>     initialize: function () 
#>     is_response_transformed: FALSE
#>     isotope_label: isotopeLabel
#>     min_peptides_protein: 2
#>     norm_value: NULL
#>     nr_children: nr_children
#>     opt_mz: 
#>     opt_rt: 
#>     opt_se: 
#>     pattern_contaminants: NULL
#>     pattern_decoys: NULL
#>     pop_response: function () 
#>     sample_name: sampleName
#>     sep: ~
#>     set_response: function (col_name) 
#>     value_vars: function () 
#>     work_intensity: abundance
#> 
if(require("SummarizedExperiment")){
   tmp <- LFQDataToSummarizedExperiment(lfqdata)
}
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
```
