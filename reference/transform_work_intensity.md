# Transform intensity column by applying a function

Transform intensity column by applying a function

## Usage

``` r
transform_work_intensity(pdata, response, .func, .funcname = NULL)
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

## Value

list with \`data\` (data.frame) and \`colname\` (new column name)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
config <- istar$config
res <- transform_work_intensity(
  istar$data, config$get_response(), .func = log2
)
#> Column added : log2_abundance
stopifnot("log2_abundance" %in% colnames(res$data))
res <- transform_work_intensity(
  istar$data, config$get_response(), .func = asinh
)
#> Column added : asinh_abundance
stopifnot("asinh_abundance" %in% colnames(res$data))
```
