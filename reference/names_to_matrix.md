# splits names and creates a matrix

splits names and creates a matrix

## Usage

``` r
names_to_matrix(names, split = "\\||\\_")
```

## Arguments

- names:

  vector with names

- split:

  patter to use to split

## Value

matrix

## See also

Other utilities:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`get_uniprot_id_from_fasta_header()`](https://wolski.github.io/prolfqua/reference/get_uniprot_id_from_fasta_header.md),
[`matrix_to_tibble()`](https://wolski.github.io/prolfqua/reference/matrix_to_tibble.md),
[`multigroup_volcano()`](https://wolski.github.io/prolfqua/reference/multigroup_volcano.md),
[`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md),
[`panel_cor()`](https://wolski.github.io/prolfqua/reference/panel_cor.md),
[`remove_na_rows()`](https://wolski.github.io/prolfqua/reference/remove_NA_rows.md),
[`table_facade()`](https://wolski.github.io/prolfqua/reference/table_facade.md)

## Examples

``` r
dat = c("bla_ra0/2_run0","bla_ra1/2_run0","bla_ra2/2_run0")
names_to_matrix(dat,split="\\_|\\/")
#>      [,1]  [,2]  [,3] [,4]  
#> [1,] "bla" "ra0" "2"  "run0"
#> [2,] "bla" "ra1" "2"  "run0"
#> [3,] "bla" "ra2" "2"  "run0"
```
