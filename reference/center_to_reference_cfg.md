# center to reference

takes the median of the lfqdareference per sample and subtracts it from
a copy of lfqdata

## Usage

``` r
center_to_reference_cfg(lfqdata, lfqdareference)
```

## Arguments

- lfqdata:

  LFQData object containing the data to center

- lfqdareference:

  LFQData object containing the reference subset

## Value

The computed result.

## Examples

``` r
# example code

bb <- sim_lfq_data_peptide_config(Nprot = 100)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
x <- LFQData$new(bb$data, bb$config)
xc <- x$get_copy()
xc$set_data(xc$data_long() |> dplyr::filter(protein_Id == "0EfVhX~3967"))
xxd <- center_to_reference_cfg(x, xc)
xxd$response()
#> [1] "centered_abundance_by_median"
x$response()
#> [1] "abundance"
```
