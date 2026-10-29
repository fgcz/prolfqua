# Table of distinct factors (sample annotation)

Table of distinct factors (sample annotation)

## Usage

``` r
table_factors(pdata, file_name, sample_name, factor_keys)
```

## Arguments

- pdata:

  data.frame

- file_name:

  character — file name column

- sample_name:

  character — sample name column

- factor_keys:

  character vector — factor column names

## See also

Other configuration:
[`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://fgcz.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://fgcz.github.io/prolfqua/reference/complete_cases.md),
[`list_to_ContrastConfiguration()`](https://fgcz.github.io/prolfqua/reference/list_to_ContrastConfiguration.md),
[`make_interaction_column()`](https://fgcz.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://fgcz.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://fgcz.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://fgcz.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors_size()`](https://fgcz.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(istar$data, istar$config)
xx <- table_factors(lfq$data_long(), lfq$file_name(), lfq$sample_name(), lfq$factor_keys())
xt <- xx |> dplyr::group_by(!!!rlang::syms(lfq$factor_keys())) |>
 dplyr::summarize(n = dplyr::n())
stopifnot(all(xt$n == 4))
```
