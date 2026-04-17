# center to reference

takes the mean or median of the lfqdareference per sample and subtracts
from lfqdata

## Usage

``` r
center_to_reference_cfg(
  lfqdata,
  lfqdareference,
  summary = c("median", "mean"),
  copy = TRUE
)
```

## Arguments

- lfqdata:

  LFQData object containing the data to center

- lfqdareference:

  LFQData object containing the reference subset

- summary:

  character, summary statistic to use ("median" or "mean")

- copy:

  logical, if TRUE return a copy, otherwise modify in place

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
xxd <- center_to_reference_cfg(x, xc, summary="median")
xxd$response()
#> [1] "centered_abundance_by_median"
xxd$data
#> NULL
center_to_reference_cfg(x, xc, summary="median", copy=FALSE)
x$response()
#> [1] "centered_abundance_by_median"
```
