# ranks precursor - peptide by intensity.

ranks precursor - peptide by intensity.

## Usage

``` r
rank_peptide_by_intensity(pdata, config)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

## Value

data.frame

## Examples

``` r


bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
res <- rank_peptide_by_intensity(bb$data, bb$config)
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank
X <-res |> dplyr::select(c(bb$config$hierarchy_keys(),
 srm_meanInt, srm_meanIntRank)) |> dplyr::distinct()
X |> dplyr::arrange(!!!rlang::syms(c(bb$config$hierarchy_keys()[1], "srm_meanIntRank"  )))
#> # A tibble: 28 × 4
#>    protein_Id  peptide_Id srm_meanInt srm_meanIntRank
#>    <chr>       <chr>            <dbl>           <int>
#>  1 0EfVhX~0087 ahQLlQY7          24.9               1
#>  2 0EfVhX~0087 ITLb4x1q          23.2               2
#>  3 0EfVhX~0087 dJkdz7so          20.7               3
#>  4 7cbcrd~5725 D5dQ4nKk          23.5               1
#>  5 9VUkAq~4703 eIC06D7g          21.1               1
#>  6 BEJI92~5282 qQ1GK8Un          23.2               1
#>  7 BEJI92~5282 HBkZvdhT          18.0               2
#>  8 CGzoYe~2147 mjHSHhoe          28.3               1
#>  9 DoWup2~5896 KVUnZ6oZ          20.5               1
#> 10 Fl4JiV~8625 wajUl0YS          25.9               1
#> # ℹ 18 more rows
```
