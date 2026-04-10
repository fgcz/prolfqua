# Table of distinct factors with group sizes

Table of distinct factors with group sizes

## Usage

``` r
table_factors_size(
  pdata,
  file_name,
  sample_name,
  factor_keys,
  factor_keys_depth
)
```

## Arguments

- pdata:

  data.frame

- file_name:

  character — file name column

- sample_name:

  character — sample name column

- factor_keys:

  character vector — all factor column names

- factor_keys_depth:

  character vector — factor columns at current depth

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
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(istar$data, istar$config)
xx <- table_factors_size(lfq$get_data(), lfq$file_name(),
  lfq$sample_name(), lfq$factor_keys(), lfq$relevant_factor_keys())
stopifnot(all(xx$n == 4))
```
