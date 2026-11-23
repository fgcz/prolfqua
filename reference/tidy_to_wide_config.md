# Transform long to wide format

Transform long to wide format

## Usage

``` r
tidy_to_wide_config(
  lfqdata,
  as.matrix = FALSE,
  file_name = FALSE,
  sep = "~lfq~",
  value = lfqdata$response(),
  row_keys = lfqdata$hierarchy_keys()
)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- as.matrix:

  if TRUE return matrix, otherwise data.frame

- file_name:

  if TRUE use file_name as column labels, otherwise sample_name

- sep:

  separator for row IDs when as.matrix = TRUE

- value:

  column name to pivot (default: lfqdata\$response())

- row_keys:

  columns identifying a row of the wide table (default: the hierarchy
  keys). They must identify one value per sample.

## Value

list with data, rowdata, and annotation (colData)

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$data_wide()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
