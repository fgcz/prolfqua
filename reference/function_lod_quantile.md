# get smallest values per sample

get smallest values per sample

## Usage

``` r
function_lod_quantile(data_matrix, percent = 10)
```

## Arguments

- data_matrix:

  numeric matrix or data.frame of abundance values

- percent:

  numeric, percentile cutoff for smallest values

## Examples

``` r
# example code

istar <- sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
xx <- lfqdata$to_wide(as.matrix=TRUE)
s <- function_lod_quantile(xx$data)
sapply(s, median)
#>     A_V1     A_V2     A_V3     A_V4     B_V1     B_V2     B_V3     B_V4 
#> 18.12727 17.26621 14.68378 18.77994 18.40038 18.48597 18.38984 17.66791 
#>  Ctrl_V1  Ctrl_V2  Ctrl_V3  Ctrl_V4 
#> 17.45284 16.73972 18.33940 17.26318 
sapply(s, mean)
#>     A_V1     A_V2     A_V3     A_V4     B_V1     B_V2     B_V3     B_V4 
#> 17.86779 16.95259 14.90349 18.53198 17.90055 18.51540 18.16227 17.62211 
#>  Ctrl_V1  Ctrl_V2  Ctrl_V3  Ctrl_V4 
#> 17.44378 16.74634 18.56473 16.92604 
```
