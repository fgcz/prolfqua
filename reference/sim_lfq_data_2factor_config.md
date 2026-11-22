# Simulate data, protein, with config with 2 factors Treatment and Background

Simulate data, protein, with config with 2 factors Treatment and
Background

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

  add missing values, default TRUE

- weight_missing:

  one number in \[0, 1\]; greater weight means more missingness

- peptide:

  if TRUE, simulate peptide-level data; if FALSE, protein-level

- seed:

  seed for reproducibility, if NULL no seed is set.

- two_factor:

  if TRUE, model `Treatment` and `Background` as two factors; if FALSE,
  as one combined factor

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
#> [1] 160   8
stopifnot("data.frame" %in% class(x$data))
stopifnot("AnalysisConfiguration" %in% class(x$config))
x <- sim_lfq_data_2factor_config(peptide = TRUE)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 7 unconfigured columns: average_prot_abundance, sd, group, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 40 rows with NA response for feature and sample combinations absent from the input.

head(x$data)
#> # A tibble: 6 × 9
#>   sample sampleName Treatment Background protein_Id  peptide_Id abundance qValue
#>   <chr>  <chr>      <chr>     <chr>      <chr>       <chr>          <dbl>  <dbl>
#> 1 A_V1   A_V1       A         X          0EfVhX~0087 ITLb4x1q        15.8      0
#> 2 A_V1   A_V1       A         X          0EfVhX~0087 ahQLlQY7        20.0      0
#> 3 A_V1   A_V1       A         X          0EfVhX~0087 dJkdz7so        15.8      0
#> 4 A_V1   A_V1       A         X          7cbcrd~5725 D5dQ4nKk        29.3      0
#> 5 A_V1   A_V1       A         X          9VUkAq~4703 eIC06D7g        20.0      0
#> 6 A_V1   A_V1       A         X          BEJI92~5282 HBkZvdhT        26.5      0
#> # ℹ 1 more variable: nrPeptides <dbl>
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
