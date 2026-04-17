# esitmate lod

esitmate lod

## Usage

``` r
estimate_lod_global(data_matrix, prop_na = 90)
```

## Arguments

- data_matrix:

  numeric matrix of abundance values

- prop_na:

  numeric, percentage threshold for NA proportion per row

## Examples

``` r
# example code
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
xx <- lfqdata$data_wide(as.matrix=TRUE)
stopifnot(length(estimate_lod_global(xx$data, prop_na = 90)) == 0)
stopifnot(length(estimate_lod_global(xx$data, prop_na = 10)) > 0)
```
