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
xc$data <- xc$data |> dplyr::filter(protein_Id == "0EfVhX~3967")
xxd <- center_to_reference_cfg(x, xc, summary="median")
xxd$response()
#> [1] "centered_abundance_by_median"
xxd$data
#> # A tibble: 4,200 × 11
#>    sample sampleName group_ isotopeLabel protein_Id  peptide_Id abundance qValue
#>    <chr>  <chr>      <chr>  <chr>        <chr>       <chr>          <dbl>  <dbl>
#>  1 A_V1   A_V1       A      light        0EfVhX~3967 IIhYJDAe        24.2      0
#>  2 A_V1   A_V1       A      light        0EfVhX~3967 SWkbauTR        25.5      0
#>  3 A_V1   A_V1       A      light        0YSKpy~2865 CBF8v3h7        16.4      0
#>  4 A_V1   A_V1       A      light        0YSKpy~2865 PdzsyJTo        NA       NA
#>  5 A_V1   A_V1       A      light        0m5WN4~6025 7uKIY8WX        12.8      0
#>  6 A_V1   A_V1       A      light        0m5WN4~6025 7xDNA2B6        19.4      0
#>  7 A_V1   A_V1       A      light        0m5WN4~6025 KT0ROM7b        22.2      0
#>  8 A_V1   A_V1       A      light        0m5WN4~6025 LYLauRlr        20.4      0
#>  9 A_V1   A_V1       A      light        0m5WN4~6025 PZ6aqY4E        16.6      0
#> 10 A_V1   A_V1       A      light        0m5WN4~6025 VaCySyEM        23.6      0
#> # ℹ 4,190 more rows
#> # ℹ 3 more variables: nr_children <dbl>, centered_abundance_by_mean <dbl>,
#> #   centered_abundance_by_median <dbl>
center_to_reference_cfg(x, xc, summary="median", copy=FALSE)
x$response()
#> [1] "centered_abundance_by_median"
```
