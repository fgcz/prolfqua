# Rank precursors or peptides by intensity

Rank precursors or peptides by intensity

## Usage

``` r
rank_peptide_by_intensity(pdata, response, hierarchy_keys)
```

## Arguments

- pdata:

  data.frame

- response:

  character - intensity column name

- hierarchy_keys:

  character vector - all hierarchy column names

## Value

data.frame

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Aggregator("topN")$aggregate()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
