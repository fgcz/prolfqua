# Simulate a two-factor design with its configuration

Simulated protein (or peptide) intensities for a 2 x 2 design of
`Treatment` and `Background` with four samples per cell, returned with a
matching
[AnalysisConfiguration](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md),
ready for `LFQData$new()`.

## Usage

``` r
sim_lfq_data_2factor_config(
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  peptide = FALSE,
  seed = 1234,
  two_factor = TRUE,
  Nprot,
  PEPTIDE,
  TWO
)
```

## Arguments

- n_proteins:

  one positive whole number of proteins to simulate

- with_missing:

  one TRUE/FALSE: add missing values (default TRUE)

- weight_missing:

  one number in \[0, 1\]; greater weight means more missingness

- peptide:

  one TRUE/FALSE: TRUE simulates peptide-level data, FALSE protein-level
  data

- seed:

  NULL or one whole number; NULL sets no seed

- two_factor:

  one TRUE/FALSE: TRUE models `Treatment` and `Background` as two
  factors, FALSE as one combined factor

- Nprot, PEPTIDE, TWO:

  deprecated; use `n_proteins`, `peptide` and `two_factor`

## Value

A list containing simulated data and an analysis configuration.

## Examples

``` r
x <- sim_lfq_data_2factor_config(peptide = FALSE)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, group, mean, N, Replicate.
#> setup_analysis: added 22 rows with NA response for feature and sample combinations absent from the input.
dim(x$data)
#> [1] 160  10
stopifnot("data.frame" %in% class(x$data))
stopifnot("AnalysisConfiguration" %in% class(x$config))
x <- sim_lfq_data_2factor_config(peptide = TRUE)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 7 unconfigured columns: average_prot_abundance, sd, group, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 40 rows with NA response for feature and sample combinations absent from the input.

head(x$data)
#> # A tibble: 6 × 11
#>   sample sampleName Treatment Background protein_Id  peptide_Id abundance
#>   <chr>  <chr>      <chr>     <chr>      <chr>       <chr>          <dbl>
#> 1 A_V1   A_V1       A         X          0EfVhX~0087 ITLb4x1q        15.8
#> 2 A_V1   A_V1       A         X          0EfVhX~0087 ahQLlQY7        20.0
#> 3 A_V1   A_V1       A         X          0EfVhX~0087 dJkdz7so        15.8
#> 4 A_V1   A_V1       A         X          7cbcrd~5725 D5dQ4nKk        29.3
#> 5 A_V1   A_V1       A         X          9VUkAq~4703 eIC06D7g        20.0
#> 6 A_V1   A_V1       A         X          BEJI92~5282 HBkZvdhT        26.5
#> # ℹ 4 more variables: identification_confidence <dbl>,
#> #   identification_confidence_feature <dbl>, nrPeptides <dbl>,
#> #   nr_children_feature <dbl>
x <- sim_lfq_data_2factor_config(peptide = TRUE, two_factor = TRUE)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 7 unconfigured columns: average_prot_abundance, sd, group, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 40 rows with NA response for feature and sample combinations absent from the input.
names(x)
#> [1] "data"   "config"
nrow(x$data) > 10
#> [1] TRUE
x$data$Treatment |> table()
#> 
#>   A   B 
#> 224 224 
```
