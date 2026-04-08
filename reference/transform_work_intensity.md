# Transform intensity

Transform intensity

## Usage

``` r
transform_work_intensity(
  pdata,
  config,
  .func,
  .funcname = NULL,
  intensity_new_name = NULL,
  deep = FALSE
)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

- .func:

  function to transform intensities e.g. log2

- .funcname:

  generates new name from name of transformation and old working
  intensity column name.

- intensity_new_name:

  column name for new intensity, default NULL

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
x <- transform_work_intensity(analysis, config, .func = log2)
#> Column added : log2_FG.Quantity
stopifnot("log2_FG.Quantity" %in% colnames(x))
config <- dd$config_f()
x <- transform_work_intensity(analysis, config, .func = asinh)
#> Column added : asinh_FG.Quantity
stopifnot("asinh_FG.Quantity" %in% colnames(x))
```
