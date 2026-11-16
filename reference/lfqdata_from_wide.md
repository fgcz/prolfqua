# Build LFQData from wide matrices

The inverse of \[tidy_to_wide_config()\]: feature x sample matrices,
their row annotation and their sample annotation back to long-format
LFQData. Row \`i\` of every matrix belongs to row \`i\` of \`rowdata\`;
matrix column names are the sample names.

## Usage

``` r
lfqdata_from_wide(values, rowdata, annotation, config, prefix = "ms")
```

## Arguments

- values:

  a matrix (features in rows, samples in columns), stored under
  \`config\$get_response()\`, or a named list of such matrices with
  identical dimensions and column names; list names are the value
  columns

- rowdata:

  data.frame with one row per matrix row and the hierarchy key columns;
  the isotope label column is optional and set to \`"light"\` if absent

- annotation:

  data.frame with one row per sample: \`config\$sample_name\` plus the
  file name and factor columns

- config:

  AnalysisConfiguration describing the result

- prefix:

  prefix passed to \[LFQData\]

## Value

LFQData

## See also

Other LFQData:
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md),
[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md),
[`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://fgcz.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://fgcz.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
lfqdata <- prolfqua::LFQData$new(dd$data, dd$config)
wide <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)
back <- lfqdata_from_wide(wide$data, wide$rowdata, wide$annotation, lfqdata$get_config())
stopifnot(nrow(back$data_long()) == nrow(lfqdata$data_long()))
```
