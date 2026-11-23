# Convert named group-mean contrasts to coefficient weights

Convert named group-mean contrasts to coefficient weights

## Usage

``` r
linfct_matrix_contrasts(linfct, contrasts, p.message = FALSE)
```

## Arguments

- linfct:

  linear functions as created by linfct_from_model

- contrasts:

  named character vector of contrasts to determine linear functions for

- p.message:

  print messages default FALSE

## Value

a matrix with one row per computable contrast and one column per model
coefficient

## Details

Rows of `linfct` represent group or interaction means and columns
represent fitted coefficients. Contrast expressions combine these rows,
so `"group_A - group_Ctrl"` becomes their row difference. Quote
interaction names containing colons with backticks. Expressions may
reference earlier named contrasts. A row naming a group unavailable in a
particular model is dropped with a warning. The returned matrix has one
row per estimable contrast and coefficient columns in model order.

## Note

Deprecated in prolfqua 1.9.5; use `Contrasts$get_linfct()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
