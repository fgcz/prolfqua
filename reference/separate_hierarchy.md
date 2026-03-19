# Separates hierarchy columns into starting columns

Separates hierarchy columns into starting columns

## Usage

``` r
separate_hierarchy(data, config)
```

## Arguments

- data:

  data.frame

- config:

  AnlalysisConfiguration

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`concrete_AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/concrete_AnalysisConfiguration.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_factors()`](https://wolski.github.io/prolfqua/reference/separate_factors.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`spread_response_by_IsotopeLabel()`](https://wolski.github.io/prolfqua/reference/spread_response_by_IsotopeLabel.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r


bb <- sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
dt <- separate_hierarchy(bb$data, bb$config)
base::setdiff(colnames(dt) ,colnames(bb$data))
#> [1] "proteinID" "idtype2"  
stopifnot(ncol(dt) >= ncol(bb$data))
```
