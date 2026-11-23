# Histogram summarizing missingness

Histogram summarizing missingness

## Usage

``` r
missingness_histogram(
  lfqdata,
  showempty = FALSE,
  factors = lfqdata$relevant_factor_keys(),
  alpha = 0.1
)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- showempty:

  show empty values

- factors:

  factor columns to use

- alpha:

  transparency

## Value

a ggplot object with density plots of the mean intensity per number of
missing values, faceted by factor levels

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Plotter()$missingness_histogram()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
