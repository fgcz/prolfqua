# Spreads isotope label heavy and light into two columns

Spreads isotope label heavy and light into two columns

## Usage

``` r
spread_response_by_IsotopeLabel(resData, config)
```

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
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r

bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
configur <- bb$config$clone(deep=TRUE)
data <- bb$data

x<-spread_response_by_IsotopeLabel(data,configur)

bb <- prolfqua_data('data_skylineSRM_HL_A')
configur <- bb$config_f()
data <- bb$analysis(bb$data, configur)
#> creating sampleName from fileName column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done

bb <- prolfqua_data('data_skylineSRM_HL_A')
conf <- bb$config_f()
analysis <- bb$analysis(bb$data, bb$config_f())
#> creating sampleName from fileName column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done
x <- spread_response_by_IsotopeLabel(analysis, conf)
```
