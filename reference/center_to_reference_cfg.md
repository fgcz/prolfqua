# Center data to a reference subset

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

invisibly, a copy of `lfqdata` whose response
`centered_abundance_by_median` is centred per sample to the median of
the reference subset.

## Examples

``` r
# example code

bb <- sim_lfq_data_peptide_config(n_proteins = 100)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 528 rows with NA response for feature and sample combinations absent from the input.
x <- LFQData$new(bb$data, bb$config)
xc <- x$get_copy()
xc$set_data(xc$data_long() |> dplyr::filter(protein_Id == "0EfVhX~3967"))
xxd <- center_to_reference_cfg(x, xc)
xxd$response()
#> [1] "centered_abundance_by_median"
x$response()
#> [1] "abundance"
```
