# Compute fold changes given Contrasts

Compute fold changes given Contrasts

## Usage

``` r
get_contrast(data, hierarchy_keys, contrasts)
```

## Arguments

- data:

  data.frame

- hierarchy_keys:

  hierarchy key columns of the AnalysisConfiguration

- contrasts:

  list of contrasts to compute

## Value

a tibble with one row per hierarchy element and contrast and columns of
the `hierarchy_keys`, `group_1`, `group_2`, `estimate`, `group_1_name`,
`group_2_name` and `contrast`

## Note

Deprecated in prolfqua 1.9.5 and will be removed; only
prolfquatakeoutmissing used it. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
