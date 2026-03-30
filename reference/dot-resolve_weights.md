# Resolve strategy weights to a matrix or vector for limma::lmFit

Shared weight resolution logic for `build_model_limma*` functions.
Handles character (column name in annotation or LFQData), matrix, or
NULL.

## Usage

``` r
.resolve_weights(lfqdata, strategy, annotation)
```

## Arguments

- lfqdata:

  an [`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md)
  object

- strategy:

  a strategy object with a `weights` field

- annotation:

  the sample-level annotation data.frame

## Value

a weight matrix, vector, or NULL
