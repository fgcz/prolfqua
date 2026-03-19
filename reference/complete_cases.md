# Complete cases

The tidy table does not need to contain missing data. This function
re-establishes the missing observations in a sample.

## Usage

``` r
complete_cases(pdata, config)
```

## Arguments

- pdata:

  data.frame

- config:

  AnlalysisConfiguration

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`concrete_AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/concrete_AnalysisConfiguration.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
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

bb <- sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done

xx <- complete_cases(bb$data, bb$config)
#> completing cases
stopifnot(nrow(bb$data) <= nrow(xx))
```
