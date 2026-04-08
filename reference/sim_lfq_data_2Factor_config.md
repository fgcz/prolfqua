# Simulate data, protein, with config with 2 factors Treatment and Background

Simulate data, protein, with config with 2 factors Treatment and
Background

## Usage

``` r
sim_lfq_data_2factor_config(
  Nprot = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  PEPTIDE = FALSE,
  seed = 1234,
  TWO = TRUE
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

- PEPTIDE:

  if TRUE, simulate peptide-level data; if FALSE, protein-level

- seed:

  seed for reproducibility, if NULL no seed is set.

- TWO:

  use two factors for modelling

## Examples

``` r
x <- sim_lfq_data_2factor_config(PEPTIDE= FALSE)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
dim(x$data)
#> [1] 160   9
stopifnot("data.frame" %in% class(x$data))
stopifnot("AnalysisConfiguration" %in% class(x$config))
x <- sim_lfq_data_2factor_config(PEPTIDE = TRUE)
#> Warning: Unknown or uninitialised column: `nr_peptides`.
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done

head(x$data)
#> # A tibble: 6 × 10
#>   sample sampleName Treatment Background isotopeLabel protein_Id  peptide_Id
#>   <chr>  <chr>      <chr>     <chr>      <chr>        <chr>       <chr>     
#> 1 A_V1   A_V1       A         X          light        0EfVhX~0087 ITLb4x1q  
#> 2 A_V1   A_V1       A         X          light        0EfVhX~0087 ahQLlQY7  
#> 3 A_V1   A_V1       A         X          light        0EfVhX~0087 dJkdz7so  
#> 4 A_V1   A_V1       A         X          light        7cbcrd~5725 D5dQ4nKk  
#> 5 A_V1   A_V1       A         X          light        9VUkAq~4703 eIC06D7g  
#> 6 A_V1   A_V1       A         X          light        BEJI92~5282 HBkZvdhT  
#> # ℹ 3 more variables: abundance <dbl>, qValue <dbl>, nr_peptides <dbl>
x <- sim_lfq_data_2factor_config(PEPTIDE = TRUE, TWO = TRUE)
#> Warning: Unknown or uninitialised column: `nr_peptides`.
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
names(x)
#> [1] "data"   "config"
nrow(x$data) > 10
#> [1] TRUE
x$data$Treatment |> table()
#> 
#>   A   B 
#> 224 224 
```
