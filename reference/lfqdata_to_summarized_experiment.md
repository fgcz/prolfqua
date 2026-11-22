# Convert an LFQData object to a SummarizedExperiment

For compatibility with Bioconductor

## Usage

``` r
lfqdata_to_summarized_experiment(lfqdata)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

## Value

SummarizedExperiment (bioconductor)

## See also

Other LFQData:
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md),
[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md),
[`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://fgcz.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md),
[`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md),
[`plot.LFQData()`](https://fgcz.github.io/prolfqua/reference/plot.LFQData.md)

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
data <- istar$data
lfqdata <- LFQData$new(data, istar$config)
lfqdata$data_wide()
#> $data
#> # A tibble: 28 × 14
#>    protein_Id peptide_Id  A_V1  A_V2  A_V3  A_V4  B_V1  B_V2  B_V3  B_V4 Ctrl_V1
#>    <chr>      <chr>      <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#>  1 0EfVhX~00… ITLb4x1q    19.4  NA    18.2  18.9  28.9  27.3  29.9  30.4    20.1
#>  2 0EfVhX~00… ahQLlQY7    26.9  25.2  26.1  25.2  26.1  26.2  24.5  NA      23.7
#>  3 0EfVhX~00… dJkdz7so    NA    16.1  14.7  NA    25.7  23.9  25.5  25.1    NA  
#>  4 7cbcrd~57… D5dQ4nKk    30.4  30.4  29.6  27.8  NA    21.6  20.9  21.0    18.8
#>  5 9VUkAq~47… eIC06D7g    17.1  17.3  19.0  18.8  16.4  NA    17.4  17.5    25.8
#>  6 BEJI92~52… HBkZvdhT    18.1  17.9  NA    NA    19.1  19.2  18.7  18.7    NA  
#>  7 BEJI92~52… qQ1GK8Un    21.7  24.0  24.7  24.8  28.6  25.8  26.2  24.0    16.0
#>  8 CGzoYe~21… mjHSHhoe    24.1  24.2  25.0  23.8  NA    31.3  31.1  31.1    28.4
#>  9 DoWup2~58… KVUnZ6oZ    24.2  24.4  23.3  23.4  18.4  18.5  18.4  17.7    19.2
#> 10 Fl4JiV~86… GsUIOl6Q    20.2  20.1  20.3  NA    25.0  25.5  24.9  25.1    17.5
#> # ℹ 18 more rows
#> # ℹ 3 more variables: Ctrl_V2 <dbl>, Ctrl_V3 <dbl>, Ctrl_V4 <dbl>
#> 
#> $annotation
#> # A tibble: 12 × 3
#>    sampleName sample  group_
#>    <chr>      <chr>   <chr> 
#>  1 A_V1       A_V1    A     
#>  2 A_V2       A_V2    A     
#>  3 A_V3       A_V3    A     
#>  4 A_V4       A_V4    A     
#>  5 B_V1       B_V1    B     
#>  6 B_V2       B_V2    B     
#>  7 B_V3       B_V3    B     
#>  8 B_V4       B_V4    B     
#>  9 Ctrl_V1    Ctrl_V1 Ctrl  
#> 10 Ctrl_V2    Ctrl_V2 Ctrl  
#> 11 Ctrl_V3    Ctrl_V3 Ctrl  
#> 12 Ctrl_V4    Ctrl_V4 Ctrl  
#> 
#> $rowdata
#> # A tibble: 28 × 2
#>    protein_Id  peptide_Id
#>    <chr>       <chr>     
#>  1 0EfVhX~0087 ITLb4x1q  
#>  2 0EfVhX~0087 ahQLlQY7  
#>  3 0EfVhX~0087 dJkdz7so  
#>  4 7cbcrd~5725 D5dQ4nKk  
#>  5 9VUkAq~4703 eIC06D7g  
#>  6 BEJI92~5282 HBkZvdhT  
#>  7 BEJI92~5282 qQ1GK8Un  
#>  8 CGzoYe~2147 mjHSHhoe  
#>  9 DoWup2~5896 KVUnZ6oZ  
#> 10 Fl4JiV~8625 GsUIOl6Q  
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
#>     hierarchy_keys: function () 
#>     hierarchy_keys_depth: function (names = TRUE) 
#>     id_vars: function () 
#>     ident_q_value: qValue
#>     ident_score: 
#>     is_response_transformed: FALSE
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
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
   tmp <- lfqdata_to_summarized_experiment(lfqdata)
}
```
