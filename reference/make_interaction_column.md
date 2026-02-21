# create interaction column from factors

create interaction column from factors

## Usage

``` r
make_interaction_column(data, columns, sep = ".")
```

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`concrete_AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/concrete_AnalysisConfiguration.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_factors()`](https://wolski.github.io/prolfqua/reference/separate_factors.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`spread_response_by_IsotopeLabel()`](https://wolski.github.io/prolfqua/reference/spread_response_by_IsotopeLabel.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
xx <- data.frame(A = c("a","a","a"), B = c("d","d","e"))
x <- make_interaction_column(xx, c("B","A"))
x <- make_interaction_column(xx, c("A"))
bb <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- bb$config
analysis <- bb$data

config$factorDepth <- 1
make_interaction_column(analysis,
   config$factor_keys_depth())
#> # A tibble: 120 × 9
#>    sample sampleName group_ isotopeLabel protein_Id abundance qValue nr_peptides
#>    <chr>  <chr>      <chr>  <chr>        <chr>          <dbl>  <dbl>       <dbl>
#>  1 A_V1   A_V1       A      light        0EfVhX~00…      20.1      0           3
#>  2 A_V1   A_V1       A      light        7cbcrd~57…      22.0      0           1
#>  3 A_V1   A_V1       A      light        9VUkAq~47…      19.8      0           1
#>  4 A_V1   A_V1       A      light        BEJI92~52…      21.2      0           2
#>  5 A_V1   A_V1       A      light        CGzoYe~21…      29.4      0           1
#>  6 A_V1   A_V1       A      light        DoWup2~58…      NA       NA          NA
#>  7 A_V1   A_V1       A      light        Fl4JiV~86…      20.1      0           4
#>  8 A_V1   A_V1       A      light        HvIpHG~90…      21.7      0           2
#>  9 A_V1   A_V1       A      light        JcKVfU~96…      34.5      0           7
#> 10 A_V1   A_V1       A      light        SGIVBl~57…      23.9      0           6
#> # ℹ 110 more rows
#> # ℹ 1 more variable: interaction <fct>
```
