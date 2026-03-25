# Vectorized version of [`linfct_matrix_contrasts`](https://wolski.github.io/prolfqua/reference/linfct_matrix_contrasts.md)

Same semantics but uses a single `dplyr::mutate(data, !!!parsed)` call
instead of one mutate per contrast. Falls back to per-expression
evaluation on error so that granular failure reporting is preserved.

## Usage

``` r
linfct_matrix_contrasts_vectorized(linfct, contrasts, p.message = FALSE)
```

## Arguments

- linfct:

  linear functions as created by linfct_from_model

- contrasts:

  named character vector of contrasts to determine linear functions for

- p.message:

  print messages default FALSE
