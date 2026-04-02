# Table of distinct factors (sample annotation)

Table of distinct factors (sample annotation)

## Usage

``` r
table_factors(pdata, configuration)
```

## Arguments

- pdata:

  data.frame

- configuration:

  AnalysisConfiguration

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r

istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done


xx <- table_factors(istar$data,istar$config )
xx
#> # A tibble: 12 × 3
#>    sample  sampleName group_
#>    <chr>   <chr>      <chr> 
#>  1 A_V1    A_V1       A     
#>  2 A_V2    A_V2       A     
#>  3 A_V3    A_V3       A     
#>  4 A_V4    A_V4       A     
#>  5 B_V1    B_V1       B     
#>  6 B_V2    B_V2       B     
#>  7 B_V3    B_V3       B     
#>  8 B_V4    B_V4       B     
#>  9 Ctrl_V1 Ctrl_V1    Ctrl  
#> 10 Ctrl_V2 Ctrl_V2    Ctrl  
#> 11 Ctrl_V3 Ctrl_V3    Ctrl  
#> 12 Ctrl_V4 Ctrl_V4    Ctrl  
xt <- xx |> dplyr::group_by(!!!rlang::syms(istar$config$factor_keys())) |>
 dplyr::summarize(n = dplyr::n())
xt
#> # A tibble: 3 × 2
#>   group_     n
#>   <chr>  <int>
#> 1 A          4
#> 2 B          4
#> 3 Ctrl       4
stopifnot(all(xt$n == 4))
```
