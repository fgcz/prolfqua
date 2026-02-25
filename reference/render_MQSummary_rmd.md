# render MQ Summary.

render MQ Summary.

## Usage

``` r
render_MQSummary_rmd(
  pdata,
  config,
  project_conf,
  pep = TRUE,
  dest_path = ".",
  dest_file_name = "QCandSampleSize.Rmd",
  workdir = tempdir(),
  format = c("pdf", "html"),
  markdown_path = c("doc/QCandSampleSize.Rmd")
)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

- project_conf:

  list with workunit_Id project_Id order_Id

- pep:

  are these peptide or protein data

- dest_path:

  destination path

- dest_file_name:

  name of pdf file

- workdir:

  working directory

- format:

  either pdf or html

- markdown_path:

  path to the Rmd template file

## See also

Other vignetteHelpers:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md)

## Examples

``` r
bb <- prolfqua_data('data_skylinePRMSample_A')
config <- bb$config_f()
analysis <- bb$analysis(bb$data, bb$config_f())
#> creating sampleName from fileName column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done
projectConfig <- list(workunit_Id = "xx", project_Id = "xy", order_Id = "z")
if(FALSE){
render_MQSummary_rmd(analysis,
  config ,
  projectConfig,
  workdir= tempdir()) # tempdir(check = FALSE))
}
```
