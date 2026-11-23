# Simulate peptide-level data with its configuration

Simulated peptide intensities of proteins in three groups (A, B, Ctrl),
with missing values that are more likely at low intensity (see
[`which_missing()`](https://fgcz.github.io/prolfqua/reference/which_missing.md)),
returned with a matching
[AnalysisConfiguration](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
ready for `LFQData$new()`.

## Usage

``` r
sim_lfq_data_peptide_config(
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  seed = 1234,
  n_per_group = 4,
  Nprot,
  N
)
```

## Arguments

- n_proteins:

  one positive whole number of proteins to simulate

- with_missing:

  one TRUE/FALSE: add missing values (default TRUE)

- weight_missing:

  one number in \[0, 1\]; greater weight means more missingness

- seed:

  NULL or one whole number; NULL sets no seed

- n_per_group:

  one positive whole number of replicates per group

- Nprot, N:

  deprecated; use `n_proteins` and `n_per_group`

## Value

A list containing simulated data and an analysis configuration.

## Examples

``` r
x <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
stopifnot("data.frame" %in% class(x$data))
stopifnot("AnalysisConfiguration" %in% class(x$config))
```
