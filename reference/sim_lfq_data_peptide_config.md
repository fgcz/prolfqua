# Simulate data, protein and peptide, with config

Simulate data, protein and peptide, with config

## Usage

``` r
sim_lfq_data_peptide_config(
  Nprot = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  seed = 1234,
  N = 4
)
```

## Arguments

- Nprot:

  one positive whole number of proteins to simulate

- with_missing:

  add missing values, default TRUE

- weight_missing:

  one number in \[0, 1\]; greater weight means more missingness

- seed:

  seed for reproducibility, if NULL no seed is set.

- N:

  number of replicates per group

## Value

A list containing simulated data and an analysis configuration.

## Examples

``` r
x <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
stopifnot("data.frame" %in% class(x$data))
stopifnot("AnalysisConfiguration" %in% class(x$config))
```
