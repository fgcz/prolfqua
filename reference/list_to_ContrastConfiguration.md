# Rebuild a ContrastConfiguration from a list

Round-trip partner of
[`R6_extract_values`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
which keeps the fields but drops the methods. Use it to restore the
column-role mapping from a serialized result artifact (e.g.
\`SummarizedExperiment\` metadata or an AnnData \`uns\` entry) so
consumers can call `has_pvalue()` and friends again.

## Usage

``` r
list_to_ContrastConfiguration(dd)
```

## Arguments

- dd:

  named list of
  [`ContrastConfiguration`](https://wolski.github.io/prolfqua/reference/ContrastConfiguration.md)
  fields

## Value

A
[`ContrastConfiguration`](https://wolski.github.io/prolfqua/reference/ContrastConfiguration.md).

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://wolski.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
cfg <- ContrastConfiguration$new(
  subject_id = "protein_Id",
  contrast_col = "Bait",
  pvalue_col = NA_character_
)
values <- R6_extract_values(cfg)
restored <- list_to_ContrastConfiguration(values)
stopifnot(all.equal(R6_extract_values(restored), values))
restored$has_pvalue()
#> [1] FALSE
```
