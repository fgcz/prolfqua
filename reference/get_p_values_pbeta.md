# Protein p-value from the p-value of the median fold-change peptide

Protein p-value from the p-value of the median fold-change peptide

## Usage

``` r
get_p_values_pbeta(median.p.value, n.obs, max.n = 10)
```

## Arguments

- median.p.value:

  median of the peptide p-values of a protein

- n.obs:

  number of peptides

- max.n:

  limit number of peptides per protein.

## Value

a numeric vector of p-values, the beta-distribution probability that the
median of `min(n.obs, max.n)` uniform p-values is at most
`median.p.value`

## Note

Deprecated in prolfqua 1.9.5; use
`build_contrast_analysis(method = "ropeca_nested")`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
