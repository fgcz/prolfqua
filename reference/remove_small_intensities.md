# Remove rows when intensity lower then threshold

Remove rows when intensity lower then threshold

## Usage

``` r
remove_small_intensities(pdata, config, threshold = 1)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

## Value

data.frame

## Examples

``` r
dd <- prolfqua_data('data_spectronautDIA250_A')
config <- dd$config_f()
analysis <- dd$analysis(dd$data,config)
#> creating sampleName from file_name column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done

config$get_response()
#> [1] "FG.Quantity"

res1 <- remove_small_intensities(analysis, config, threshold=1 )
res1000 <- remove_small_intensities(analysis, config, threshold=1000 )
stopifnot(nrow(res1) >  nrow(res1000))
```
