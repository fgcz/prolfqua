# Remove rows when intensity lower then threshold

Remove rows when intensity lower then threshold

## Usage

``` r
remove_small_intensities(pdata, response, threshold = 1)
```

## Arguments

- pdata:

  data.frame

- response:

  character — name of the intensity column

- threshold:

  numeric — minimum intensity to keep (default 1)

## Value

data.frame

## Examples

``` r
istar <- sim_lfq_data_peptide_config(Nprot = 20)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
res1 <- remove_small_intensities(lfqdata$data_long(), lfqdata$response(), threshold = 1)
res1000 <- remove_small_intensities(lfqdata$data_long(), lfqdata$response(), threshold = 1000)
stopifnot(nrow(res1) > nrow(res1000))
```
