# Summarize hierarchy counts

E.g compute number of peptides for each protein: for each entry of
`hierarchy` (and of `factors`), the number of distinct entries of the
other hierarchy levels, in columns named `<level>_n`. Counts
experiment-wide, whether observed or not.

## Usage

``` r
summarize_hierarchy(
  pdata,
  hierarchy_keys,
  hierarchy = hierarchy_keys,
  factors = character()
)
```

## Arguments

- pdata:

  data.frame

- hierarchy_keys:

  character vector - all hierarchy column names

- hierarchy:

  character vector - hierarchy level to group by (default
  hierarchy_keys)

- factors:

  character vector - factor columns to include

## Value

a tibble with one row per combination of the `factors` and `hierarchy`
columns and, for each other hierarchy level, a column `<level>_n` with
the number of its distinct entries

## Note

Deprecated in prolfqua 1.9.5; use `LFQData$summarize_hierarchy()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
