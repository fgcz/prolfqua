# transform long to wide

transform long to wide

## Usage

``` r
tidy_to_wide_config(
  data,
  config,
  as.matrix = FALSE,
  fileName = FALSE,
  sep = "~lfq~",
  value = config$get_response()
)
```

## Value

list with data, rowdata, and annotation (colData)

## Examples

``` r

dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- dd$config
data <- dd$data
res <- tidy_to_wide_config(data, config)
testthat::expect_equal(nrow(res$rowdata), nrow(res$data))
testthat::expect_equal(ncol(res$data) - ncol(res$rowdata) , nrow(res$annotation))
res <- tidy_to_wide_config(data, config, as.matrix = TRUE)
stopifnot(all(dim(res$data) == c(28,  12)))
stopifnot(all(dim(res$annotation) == c(12,  4)))
stopifnot(all(dim(res$rowdata) == c(28, 3)))

res <- scale(res$data)
tidy_to_wide_config(data, config,  value = config$nr_children)
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


xt <- prolfqua::LFQData$new(dd$data, dd$config)
xt$data$nr_children
#>   [1]  1  1 NA  1  1  1  1  1  1  1  1 NA  1  1 NA  1  1  1  1  1  1  1  1  1  1
#>  [26] NA  1 NA NA  1  1  1  1  1  1  1  1  1  1 NA  1  1  1  1  1  1  1  1  1  1
#>  [51]  1  1  1  1  1  1  1  1  1  1  1 NA  1  1  1  1  1  1  1 NA  1  1  1  1  1
#>  [76]  1  1  1  1  1  1  1  1 NA  1  1 NA  1  1 NA  1  1  1 NA  1 NA  1  1  1  1
#> [101]  1  1  1  1  1  1  1  1  1  1  1 NA  1  1  1 NA  1  1  1 NA  1  1  1  1  1
#> [126]  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 NA  1  1  1  1  1
#> [151]  1  1  1  1  1  1  1  1  1  1  1  1 NA  1  1  1  1  1  1  1  1  1  1  1  1
#> [176]  1  1  1  1 NA  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 NA  1  1
#> [201]  1  1  1  1  1  1 NA  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
#> [226]  1 NA  1  1 NA  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 NA  1  1 NA
#> [251]  1  1  1 NA  1  1  1  1 NA  1  1 NA  1  1  1  1  1  1  1  1  1  1  1  1  1
#> [276]  1  1  1  1  1 NA  1  1 NA  1 NA  1  1  1  1  1  1  1  1  1  1  1  1  1  1
#> [301]  1  1  1  1  1 NA  1  1  1  1  1  1  1  1  1  1  1  1 NA  1  1  1 NA  1  1
#> [326]  1  1  1  1  1 NA  1  1  1  1  1
#xt$config$is_response_transformed <- TRUE
res <- xt$get_Aggregator()
x <- res$medpolish()
#> Warning: You did not transform the intensities.medpolish works best with already variance stabilized intensities.Use LFQData$get_Transformer to transform the data :abundance
#> starting aggregation
towide <- tidy_to_wide_config(x$data, x$config,  value = x$config$nr_children)

dd <- prolfqua::sim_lfq_data_protein_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
dd$config$nr_children
#> [1] "nr_peptides"
dd$data
#> # A tibble: 120 × 8
#>    sample sampleName group_ isotopeLabel protein_Id abundance qValue nr_peptides
#>    <chr>  <chr>      <chr>  <chr>        <chr>          <dbl>  <dbl>       <dbl>
#>  1 A_V1   A_V1       A      light        0EfVhX~00…      20.1      0           3
#>  2 A_V1   A_V1       A      light        7cbcrd~57…      22.0      0           1
#>  3 A_V1   A_V1       A      light        9VUkAq~47…      19.8      0           1
#>  4 A_V1   A_V1       A      light        BEJI92~52…      21.2      0           2
#>  5 A_V1   A_V1       A      light        CGzoYe~21…      29.4      0           1
#>  6 A_V1   A_V1       A      light        DoWup2~58…      NA       NA          NA
#>  7 A_V1   A_V1       A      light        Fl4JiV~86…      20.1      0           4
#>  8 A_V1   A_V1       A      light        HvIpHG~90…      21.7      0           2
#>  9 A_V1   A_V1       A      light        JcKVfU~96…      34.5      0           7
#> 10 A_V1   A_V1       A      light        SGIVBl~57…      23.9      0           6
#> # ℹ 110 more rows
xt <- tidy_to_wide_config(dd$data, dd$config,  value = dd$config$nr_children)
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
```
