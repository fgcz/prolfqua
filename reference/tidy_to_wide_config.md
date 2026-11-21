# Transform long to wide format

Transform long to wide format

## Usage

``` r
tidy_to_wide_config(
  lfqdata,
  as.matrix = FALSE,
  file_name = FALSE,
  sep = "~lfq~",
  value = lfqdata$response(),
  row_keys = c(lfqdata$hierarchy_keys(), lfqdata$isotope_label())
)
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- as.matrix:

  if TRUE return matrix, otherwise data.frame

- file_name:

  if TRUE use file_name as column labels, otherwise sample_name

- sep:

  separator for row IDs when as.matrix = TRUE

- value:

  column name to pivot (default: lfqdata\$response())

- row_keys:

  columns identifying a row of the wide table (default: the hierarchy
  keys and the isotope label). They must identify one value per sample;
  e.g. drop the isotope label only for single-label data.

## Value

list with data, rowdata, and annotation (colData)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- prolfqua::LFQData$new(dd$data, dd$config)
res <- tidy_to_wide_config(lfqdata)
testthat::expect_equal(nrow(res$rowdata), nrow(res$data))
testthat::expect_equal(ncol(res$data) - ncol(res$rowdata) , nrow(res$annotation))
res <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)
stopifnot(all(dim(res$data) == c(28,  12)))
stopifnot(all(dim(res$annotation) == c(12,  4)))
stopifnot(all(dim(res$rowdata) == c(28, 3)))

res <- scale(res$data)
tidy_to_wide_config(lfqdata, value = lfqdata$nr_children_col())
#> $data
#> # A tibble: 28 × 15
#>    protein_Id  peptide_Id isotopeLabel  A_V1  A_V2  A_V3  A_V4  B_V1  B_V2  B_V3
#>    <chr>       <chr>      <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 0EfVhX~0087 ITLb4x1q   light            1    NA     1     1     1     1     1
#>  2 0EfVhX~0087 ahQLlQY7   light            1     1     1     1     1     1     1
#>  3 0EfVhX~0087 dJkdz7so   light           NA     1     1    NA     1     1     1
#>  4 7cbcrd~5725 D5dQ4nKk   light            1     1     1     1    NA     1     1
#>  5 9VUkAq~4703 eIC06D7g   light            1     1     1     1     1    NA     1
#>  6 BEJI92~5282 HBkZvdhT   light            1     1    NA    NA     1     1     1
#>  7 BEJI92~5282 qQ1GK8Un   light            1     1     1     1     1     1     1
#>  8 CGzoYe~2147 mjHSHhoe   light            1     1     1     1    NA     1     1
#>  9 DoWup2~5896 KVUnZ6oZ   light            1     1     1     1     1     1     1
#> 10 Fl4JiV~8625 GsUIOl6Q   light            1     1     1    NA     1     1     1
#> # ℹ 18 more rows
#> # ℹ 5 more variables: B_V4 <dbl>, Ctrl_V1 <dbl>, Ctrl_V2 <dbl>, Ctrl_V3 <dbl>,
#> #   Ctrl_V4 <dbl>
#> 
#> $annotation
#> # A tibble: 12 × 4
#>    sampleName sample  group_ isotopeLabel
#>    <chr>      <chr>   <chr>  <chr>       
#>  1 A_V1       A_V1    A      light       
#>  2 A_V2       A_V2    A      light       
#>  3 A_V3       A_V3    A      light       
#>  4 A_V4       A_V4    A      light       
#>  5 B_V1       B_V1    B      light       
#>  6 B_V2       B_V2    B      light       
#>  7 B_V3       B_V3    B      light       
#>  8 B_V4       B_V4    B      light       
#>  9 Ctrl_V1    Ctrl_V1 Ctrl   light       
#> 10 Ctrl_V2    Ctrl_V2 Ctrl   light       
#> 11 Ctrl_V3    Ctrl_V3 Ctrl   light       
#> 12 Ctrl_V4    Ctrl_V4 Ctrl   light       
#> 
#> $rowdata
#> # A tibble: 28 × 3
#>    protein_Id  peptide_Id isotopeLabel
#>    <chr>       <chr>      <chr>       
#>  1 0EfVhX~0087 ITLb4x1q   light       
#>  2 0EfVhX~0087 ahQLlQY7   light       
#>  3 0EfVhX~0087 dJkdz7so   light       
#>  4 7cbcrd~5725 D5dQ4nKk   light       
#>  5 9VUkAq~4703 eIC06D7g   light       
#>  6 BEJI92~5282 HBkZvdhT   light       
#>  7 BEJI92~5282 qQ1GK8Un   light       
#>  8 CGzoYe~2147 mjHSHhoe   light       
#>  9 DoWup2~5896 KVUnZ6oZ   light       
#> 10 Fl4JiV~8625 GsUIOl6Q   light       
#> # ℹ 18 more rows
#> 

res <- lfqdata$get_Aggregator("medpolish")
#> Warning: You did not transform the intensities. medpolish works best with already variance stabilized intensities. Use LFQData$get_Transformer to transform the data: abundance
x <- res$aggregate()
towide <- tidy_to_wide_config(x, value = x$nr_children_col())

dd <- prolfqua::sim_lfq_data_protein_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 5 unconfigured columns: average_prot_abundance, sd, mean, N, Replicate.
#> setup_analysis: added 16 rows with NA response for feature and sample combinations absent from the input.
lfqprot <- prolfqua::LFQData$new(dd$data, dd$config)
xt <- tidy_to_wide_config(lfqprot, value = lfqprot$nr_children_col())
xt$data
#> # A tibble: 10 × 14
#>    protein_Id  isotopeLabel  A_V1  A_V2  A_V3  A_V4  B_V1  B_V2  B_V3  B_V4
#>    <chr>       <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 0EfVhX~0087 light            3     3    NA     3     3     3     3     3
#>  2 7cbcrd~5725 light            1     1     1     1    NA    NA     1    NA
#>  3 9VUkAq~4703 light            1     1     1     1     1     1     1     1
#>  4 BEJI92~5282 light            2    NA     2     2     2     2     2     2
#>  5 CGzoYe~2147 light            1     1     1     1     1     1     1     1
#>  6 DoWup2~5896 light           NA     1     1     1     1     1    NA     1
#>  7 Fl4JiV~8625 light            4     4    NA     4     4     4     4     4
#>  8 HvIpHG~9079 light            2     2    NA     2     2     2     2     2
#>  9 JcKVfU~9653 light            7     7     7     7     7     7     7     7
#> 10 SGIVBl~5782 light            6     6     6     6     6     6    NA     6
#> # ℹ 4 more variables: Ctrl_V1 <dbl>, Ctrl_V2 <dbl>, Ctrl_V3 <dbl>,
#> #   Ctrl_V4 <dbl>

# row names without the isotope label
xh <- tidy_to_wide_config(lfqprot, as.matrix = TRUE, row_keys = lfqprot$hierarchy_keys())
stopifnot(!any(grepl("light", rownames(xh$data))))
```
