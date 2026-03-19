# get the difference of two dataset where one is a subset of the other.

get the difference of two dataset where one is a subset of the other.

## Usage

``` r
filter_difference(x, y, config)
```

## Arguments

- x:

  data.frame

- y:

  data.frame

- config:

  AnlysisConfiguration

## Value

data.frame

## Examples

``` r


istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
istar$config <- istar$config
istar_data <- istar$data
filterPep <- prolfqua:::filter_proteins_by_peptide_count( istar_data ,  istar$config )
#> Column added : nr_peptide_Id_IN_protein_Id
tmp <- filter_difference(istar_data, filterPep$data, istar$config)
stopifnot(nrow(istar_data )  - nrow(filterPep$data) == nrow(tmp))
tmp <- filter_difference(filterPep$data, istar_data , istar$config)
stopifnot(nrow(istar_data )  - nrow(filterPep$data) == nrow(tmp))
```
