# normal pairs plot with different pch and plus abline

normal pairs plot with different pch and plus abline

## Usage

``` r
pairs_w_abline(dataframe, legend = FALSE, pch = ".", ...)
```

## Arguments

- dataframe:

  data matrix or data.frame as normally passed to pairs

- legend:

  add legend to plots

- pch:

  point type default "."

- ...:

  params usually passed to pairs

## See also

also [`pairs`](https://rdrr.io/r/graphics/pairs.html)

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_UniprotID_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_UniprotID_from_fasta_header.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`panel_hist()`](https://wolski.github.io/prolfqua/reference/panel_hist.md),
[`remove_NA_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
tmp = matrix(rep((1:100),times = 4) + rnorm(100*4,0,3),ncol=4)
pairs_w_abline(tmp,log="xy",main="small data")
#> Warning: 2 y values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 2 x values <= 0 omitted from logarithmic plot
#> Warning: 2 x values <= 0 omitted from logarithmic plot
#> Warning: 2 y values <= 0 omitted from logarithmic plot
#> Warning: 2 x values <= 0 omitted from logarithmic plot
#> Warning: 2 x values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 2 y values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 2 y values <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot

pairs_w_abline(tmp,log="xy",main="small data", legend=TRUE)
#> Warning: 2 y values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 2 x values <= 0 omitted from logarithmic plot
#> Warning: 2 x values <= 0 omitted from logarithmic plot
#> Warning: 2 y values <= 0 omitted from logarithmic plot
#> Warning: 2 x values <= 0 omitted from logarithmic plot
#> Warning: 2 x values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 2 y values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 2 y values <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
```
