# Rebuild a ContrastConfiguration from a list

Round-trip partner of
[`r6_extract_values`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md),
which keeps the fields but drops the methods. Use it to restore the
column-role mapping from a serialized result artifact (e.g.
`SummarizedExperiment` metadata or an AnnData `uns` entry) so consumers
can call `has_pvalue()` and friends again.

## Usage

``` r
list_to_contrast_configuration(dd)
```

## Arguments

- dd:

  named list of
  [`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md)
  fields

## Value

A
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md).

## See also

Other configuration:
[`AnalysisConfiguration`](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`ContrastConfiguration`](https://fgcz.github.io/prolfqua/reference/ContrastConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`r6_extract_values()`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md),
[`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md)

## Examples

``` r
cfg <- ContrastConfiguration$new(
  subject_id = "protein_Id",
  contrast_col = "Bait",
  pvalue_col = NA_character_
)
values <- r6_extract_values(cfg)
restored <- list_to_contrast_configuration(values)
stopifnot(all.equal(r6_extract_values(restored), values))
restored$has_pvalue()
#> [1] FALSE
```
