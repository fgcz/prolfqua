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

  number of proteins to simulate

- with_missing:

  add missing values, default TRUE

- weight_missing:

  controls proportion of missing values; greater weight means more
  missingness

- seed:

  seed for reproducibility, if NULL no seed is set.

- paired:

  if TRUE, add a paired subject factor to the design

## Examples

``` r
x <- sim_lfq_data_protein_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
stopifnot("data.frame" %in% class(x$data))
stopifnot("AnalysisConfiguration" %in% class(x$config))
x <- sim_lfq_data_protein_config(with_missing = FALSE)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done

stopifnot(sum(is.na(x$data$abundance)) == 0)
xp <- sim_lfq_data_protein_config(with_missing = FALSE, paired = TRUE)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
stopifnot(length(xp$config$factors) == 2)
stopifnot(nrow(xp$data) == nrow(x$data))
```
