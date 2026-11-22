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
  row_keys = lfqdata$hierarchy_keys()
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
  keys). They must identify one value per sample.

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
stopifnot(nrow(res$rowdata) == nrow(res$data))
stopifnot(ncol(res$data) - ncol(res$rowdata) == nrow(res$annotation))
res <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)
stopifnot(all(dim(res$data) == c(28,  12)))
stopifnot(all(dim(res$annotation) == c(12,  3)))
stopifnot(all(dim(res$rowdata) == c(28, 2)))

res <- scale(res$data)
tidy_to_wide_config(lfqdata, value = lfqdata$nr_children_col())
#> $data
#> # A tibble: 28 × 14
#>    protein_Id peptide_Id  A_V1  A_V2  A_V3  A_V4  B_V1  B_V2  B_V3  B_V4 Ctrl_V1
#>    <chr>      <chr>      <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>
#>  1 0EfVhX~00… ITLb4x1q       1    NA     1     1     1     1     1     1       1
#>  2 0EfVhX~00… ahQLlQY7       1     1     1     1     1     1     1    NA       1
#>  3 0EfVhX~00… dJkdz7so      NA     1     1    NA     1     1     1     1      NA
#>  4 7cbcrd~57… D5dQ4nKk       1     1     1     1    NA     1     1     1       1
#>  5 9VUkAq~47… eIC06D7g       1     1     1     1     1    NA     1     1       1
#>  6 BEJI92~52… HBkZvdhT       1     1    NA    NA     1     1     1     1      NA
#>  7 BEJI92~52… qQ1GK8Un       1     1     1     1     1     1     1     1       1
#>  8 CGzoYe~21… mjHSHhoe       1     1     1     1    NA     1     1     1       1
#>  9 DoWup2~58… KVUnZ6oZ       1     1     1     1     1     1     1     1       1
#> 10 Fl4JiV~86… GsUIOl6Q       1     1     1    NA     1     1     1     1       1
#> # ℹ 18 more rows
#> # ℹ 3 more variables: Ctrl_V2 <dbl>, Ctrl_V3 <dbl>, Ctrl_V4 <dbl>
#> 
#> $annotation
#> # A tibble: 12 × 3
#>    sampleName sample  group_
#>    <chr>      <chr>   <chr> 
#>  1 A_V1       A_V1    A     
#>  2 A_V2       A_V2    A     
#>  3 A_V3       A_V3    A     
#>  4 A_V4       A_V4    A     
#>  5 B_V1       B_V1    B     
#>  6 B_V2       B_V2    B     
#>  7 B_V3       B_V3    B     
#>  8 B_V4       B_V4    B     
#>  9 Ctrl_V1    Ctrl_V1 Ctrl  
#> 10 Ctrl_V2    Ctrl_V2 Ctrl  
#> 11 Ctrl_V3    Ctrl_V3 Ctrl  
#> 12 Ctrl_V4    Ctrl_V4 Ctrl  
#> 
#> $rowdata
#> # A tibble: 28 × 2
#>    protein_Id  peptide_Id
#>    <chr>       <chr>     
#>  1 0EfVhX~0087 ITLb4x1q  
#>  2 0EfVhX~0087 ahQLlQY7  
#>  3 0EfVhX~0087 dJkdz7so  
#>  4 7cbcrd~5725 D5dQ4nKk  
#>  5 9VUkAq~4703 eIC06D7g  
#>  6 BEJI92~5282 HBkZvdhT  
#>  7 BEJI92~5282 qQ1GK8Un  
#>  8 CGzoYe~2147 mjHSHhoe  
#>  9 DoWup2~5896 KVUnZ6oZ  
#> 10 Fl4JiV~8625 GsUIOl6Q  
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
#> # A tibble: 10 × 13
#>    protein_Id   A_V1  A_V2  A_V3  A_V4  B_V1  B_V2  B_V3  B_V4 Ctrl_V1 Ctrl_V2
#>    <chr>       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl>
#>  1 0EfVhX~0087     3     3    NA     3     3     3     3     3       3       3
#>  2 7cbcrd~5725     1     1     1     1    NA    NA     1    NA       1       1
#>  3 9VUkAq~4703     1     1     1     1     1     1     1     1      NA      NA
#>  4 BEJI92~5282     2    NA     2     2     2     2     2     2       2       2
#>  5 CGzoYe~2147     1     1     1     1     1     1     1     1       1       1
#>  6 DoWup2~5896    NA     1     1     1     1     1    NA     1      NA      NA
#>  7 Fl4JiV~8625     4     4    NA     4     4     4     4     4       4       4
#>  8 HvIpHG~9079     2     2    NA     2     2     2     2     2       2       2
#>  9 JcKVfU~9653     7     7     7     7     7     7     7     7       7       7
#> 10 SGIVBl~5782     6     6     6     6     6     6    NA     6       6       6
#> # ℹ 2 more variables: Ctrl_V3 <dbl>, Ctrl_V4 <dbl>
```
