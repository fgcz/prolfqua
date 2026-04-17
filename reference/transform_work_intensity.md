# Transform intensity column by applying a function

Transform intensity column by applying a function

## Usage

``` r
transform_work_intensity(
  pdata,
  response,
  .func,
  .funcname = NULL,
  intensity_new_name = NULL
)
```

## Arguments

- pdata:

  data.frame

- response:

  character, name of the response column to transform

- .func:

  function to transform intensities e.g. log2

- .funcname:

  name of function (used for creating new column name)

- intensity_new_name:

  column name for new intensity, default NULL

## Value

list with \`data\` (data.frame) and \`colname\` (new column name)

## Examples

``` r
dd <- prolfqua_data('data_spectronautDIA250_A')
config <- dd$config_f()
analysis <- dd$analysis(dd$data, config)
#> creating sampleName from file_name column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done
res <- transform_work_intensity(
  analysis, config$get_response(), .func = log2
)
#> Column added : log2_FG.Quantity
stopifnot("log2_FG.Quantity" %in% colnames(res$data))
config <- dd$config_f()
res <- transform_work_intensity(
  analysis, config$get_response(), .func = asinh
)
#> Column added : asinh_FG.Quantity
stopifnot("asinh_FG.Quantity" %in% colnames(res$data))
```
