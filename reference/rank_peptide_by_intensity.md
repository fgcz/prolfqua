# ranks precursor - peptide by intensity.

ranks precursor - peptide by intensity.

## Usage

``` r
rank_peptide_by_intensity(pdata, response, hierarchy_keys)
```

## Arguments

- pdata:

  data.frame

- response:

  character — intensity column name

- hierarchy_keys:

  character vector — all hierarchy column names

## Value

data.frame

## Examples

``` r
bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq <- LFQData$new(bb$data, bb$config)
res <- rank_peptide_by_intensity(lfq$data_long(), lfq$response(), lfq$hierarchy_keys())
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank
X <- res |> dplyr::select(c(lfq$hierarchy_keys(),
 srm_meanInt, srm_meanIntRank)) |> dplyr::distinct()
X |> dplyr::arrange(!!!rlang::syms(c(lfq$hierarchy_keys()[1], "srm_meanIntRank")))
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
