# smoothScatter pairs

smoothScatter pairs

## Usage

``` r
pairs_smooth(data, legend = FALSE, ...)
```

## Arguments

- data:

  data matrix or data.frame as normally passed to pairs

- legend:

  add legend to plots

- ...:

  params usually passed to pairs

## See also

also [`pairs`](https://rdrr.io/r/graphics/pairs.html)

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_uniprot_id_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`names_to_matrix()`](https://wolski.github.io/prolfqua/reference/names_to_matrix.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`remove_na_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
tmp = matrix(rep((1:100),times = 4) + rnorm(100*4,0,3),ncol=4)
pairs_smooth(tmp,main="small data", legend=TRUE)

pairs_smooth(tmp,main="small data")

pairs_smooth(tmp,log="xy",main="small data", legend=TRUE)
#> Warning: 3 y values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 3 x values <= 0 omitted from logarithmic plot
#> Warning: 3 x values <= 0 omitted from logarithmic plot
#> Warning: 3 y values <= 0 omitted from logarithmic plot
#> Warning: 3 x values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 3 x values <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 3 y values <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 3 y values <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
#> Warning: 1 x value <= 0 omitted from logarithmic plot
#> Warning: 1 y value <= 0 omitted from logarithmic plot
```
