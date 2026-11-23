# Adjust p-values for multiple testing

Benjamini-Hochberg false discovery rate, computed separately within each
contrast (or other group).

## Usage

``` r
adjust_p_values(
  contrast_df,
  column = "p.value",
  group_by_col = "contrast",
  newname = "FDR"
)
```

## Arguments

- contrast_df:

  data.frame with p-values to adjust

- column:

  name of column containing p-values

- group_by_col:

  column(s) to group by before adjusting (e.g. contrast), or NULL for no
  grouping

- newname:

  name of the new column with adjusted p-values

## Value

`contrast_df`, ungrouped, with the Benjamini-Hochberg adjusted p-values
in the column `newname`, adjusted within each group of `group_by_col`.

## Examples

``` r
bb <- c(runif(1000), rexp(1500,rate=5))
length(bb)
#> [1] 2500
bb <- bb[bb < 1]
length(bb)
#> [1] 2492
bb <- bb[1:2000]
hist(bb)

data <- data.frame(contrast = rep(LETTERS[1:5],400), p.value = bb)

dataX <- adjust_p_values(data)
Adata <- dataX |> dplyr::filter(contrast == "A")
stopifnot(all.equal(Adata$FDR, p.adjust(Adata$p.value, method="BH")))
data2 <- adjust_p_values(data, group_by_col = NULL)
stopifnot(all.equal(data2$FDR, p.adjust(data2$p.value, method="BH")))

```
