# Simulate data, protein, with config

Simulate data, protein, with config

## Usage

``` r
sim_lfq_data_protein_config(
  Nprot = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  seed = 1234,
  paired = FALSE
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

- paired:

  if TRUE, add a paired subject factor to the design

## Value

A list containing simulated data and an analysis configuration.

## Examples

``` r
x <- sim_lfq_data_protein_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 16 rows with NA response for feature and sample combinations absent from the input.
stopifnot("data.frame" %in% class(x$data))
stopifnot("AnalysisConfiguration" %in% class(x$config))
x <- sim_lfq_data_protein_config(with_missing = FALSE)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.

stopifnot(sum(is.na(x$data$abundance)) == 0)
xp <- sim_lfq_data_protein_config(with_missing = FALSE, paired = TRUE)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
stopifnot(length(xp$config$factors) == 2)
stopifnot(nrow(xp$data) == nrow(x$data))
```
