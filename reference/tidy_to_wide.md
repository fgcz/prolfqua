# Transform tidy table into a table with a column of responses for each sample

Transform tidy table into a table with a column of responses for each
sample

## Usage

``` r
tidy_to_wide(data, row_ids, column_labels, value)
```

## Examples

``` r
pdata <- data.frame(
  protein_Id = c("P1", "P1", "P2", "P2"),
  sampleName = c("S1", "S2", "S1", "S2"),
  abundance = c(10, 12, 20, 25)
)
tidy_to_wide(pdata, row_ids = "protein_Id", column_labels = "sampleName", value = "abundance")
#> # A tibble: 2 × 3
#>   protein_Id    S1    S2
#>   <chr>      <dbl> <dbl>
#> 1 P1            10    12
#> 2 P2            20    25
```
