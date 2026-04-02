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

  number of proteins to simulate

- with_missing:

  add missing values, default TRUE

- weight_missing:

  controls proportion of missing values; greater weight means more
  missingness

- seed:

  seed for reproducibility, if NULL no seed is set.

- N:

  number of replicates per group

## Examples

``` r
x <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
stopifnot("data.frame" %in% class(x$data))
stopifnot("AnalysisConfiguration" %in% class(x$config))
```
