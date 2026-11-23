# Create an ANOVA extractor returning a data frame

Create an ANOVA extractor returning a data frame

## Usage

``` r
get_anova_df(test = "F")
```

## Arguments

- test:

  one case-sensitive test choice: "F", "Chisq", "LRT", "Rao" or "Cp".
  The fitted model's ANOVA method must support the chosen test.

## Value

An
[`AnovaExtractor`](https://fgcz.github.io/prolfqua/reference/AnovaExtractor.md)
object.

## Note

Deprecated in prolfqua 1.9.5; use `Model$get_anova()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
