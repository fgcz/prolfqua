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

list with \`data\` (data.frame), \`colname\` (new column name) and
\`zeros\` (logical, rows whose zero intensity became \`NA\`)

## Details

For \`log\`, \`log2\` and \`log10\`, zero intensities mean "not
detected" and become \`NA\`; LFQData records them in its provenance.
Negative intensities give \`NaN\` with a warning, which fitting and
aggregation then reject.

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
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
