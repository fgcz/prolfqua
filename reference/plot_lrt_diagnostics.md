# Write LRT diagnostic plots to PDF

Writes a histogram of p-values and an empirical CDF to a PDF file.
Called by
[`LR_test`](https://wolski.github.io/prolfqua/reference/LR_test.md) when
`path` is non-NULL.

## Usage

``` r
plot_lrt_diagnostics(result, model_name, model_name_int, path)
```

## Arguments

- result:

  tibble with a `likelihood_ratio_test.pValue` column

- model_name:

  name of the full model

- model_name_int:

  name of the reduced/interaction model

- path:

  directory to write the PDF into
