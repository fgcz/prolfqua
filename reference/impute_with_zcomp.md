# Impute missing values using zCompositions

Replaces missing values in an LFQData object using methods from the
zCompositions package. The limit of detection (LOD) can be estimated
globally or per-sample using quantiles.

## Usage

``` r
impute_with_zcomp(
  lfqdata,
  method = c("multRepl", "GBM", "SQ", "BL", "CZM"),
  lod = c("global", "quantile")
)
```

## Arguments

- lfqdata:

  LFQData object containing the data to impute

- method:

  imputation method: "multRepl" (multiplicative replacement), "GBM",
  "SQ", "BL", or "CZM" (passed to zCompositions)

- lod:

  limit of detection strategy, either "global" or "quantile"

## Value

the modified LFQData object (lfqdata), with imputed values

## Note

This function assumes that missing values are Missing Completely At
Random (MCAR) or Missing At Random (MAR). If missingness is
abundance-dependent (MNAR, common in proteomics DDA), the imputed values
and downstream statistics may be biased. For MNAR-aware analysis,
consider packages such as proDA or DEqMS.

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(dd$data, dd$config)
if (requireNamespace("zCompositions", quietly = TRUE)) {
  wide_before <- lfqdata$to_wide(as.matrix = TRUE)
  has_na_before <- any(is.na(wide_before$data))
  lfqdata <- impute_with_zcomp(lfqdata, method = "multRepl", lod = "global")
  wide_after <- lfqdata$to_wide(as.matrix = TRUE)
  has_na_after <- any(is.na(wide_after$data))
  stopifnot(has_na_before || !has_na_after)
  stopifnot(!has_na_after)
}
#> Warning: Expected 3 pieces. Missing pieces filled with `NA` in 336 rows [1, 2, 3, 4, 5,
#> 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id, peptide_Id)`
```
