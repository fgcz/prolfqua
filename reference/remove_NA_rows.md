# Removes rows with more than thresh NA's from matrix

Removes rows with more than thresh NA's from matrix

## Usage

``` r
remove_NA_rows(obj, thresh = 0)
```

## Arguments

- obj:

  matrix or dataframe

- thresh:

  \- maximum number of NA's / row - if more the row will be removed

## Value

matrix

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_UniprotID_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_UniprotID_from_fasta_header.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`pairs_w_abline()`](https://wolski.github.io/prolfqua/reference/pairs_w_abline.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`panel_hist()`](https://wolski.github.io/prolfqua/reference/panel_hist.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
obj = matrix(rnorm(10*10),ncol=10)
dim(obj)
#> [1] 10 10
obj[3,3] = NA
x1 = remove_NA_rows(obj, thresh=0)
stopifnot(all(c(9,10)==dim(x1)))
x2 = remove_NA_rows(obj, thresh=1)
stopifnot(all(c(10,10)==dim(x2)))
```
