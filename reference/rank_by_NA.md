# Ranks peptides/precursors of a protein by NAs (adds new column .NARank)

Ranks peptides/precursors of a protein by NAs (adds new column .NARank)

## Usage

``` r
rank_by_NA(pdata, config)
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

bb <- prolfqua::prolfqua_data('data_spectronautDIA250_A')
config <- bb$config_f()
analysis <- bb$analysis(bb$data, bb$config_f())
#> creating sampleName from fileName column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done
res <- rank_by_NA(analysis, config)
#> Joining with `by = join_by(protein_Id, peptide_Id,
#> modPeptide_Id, precursor_Id)`
#> Columns added : srm_NrNotNAs srm_NrNotNARank
colnames(res)
#>  [1] "R.FileName"      "sampleName"      "coding"          "sex"            
#>  [5] "age"             "Sample_id"       "Isotope.Label"   "protein_Id"     
#>  [9] "peptide_Id"      "modPeptide_Id"   "precursor_Id"    "FG.Quantity"    
#> [13] "EG.Qvalue"       "nr_children"     "srm_NrNotNAs"    "srm_NrNotNARank"
x <- res |>
  dplyr::select(config$hierarchy_keys()[1],
    config$hierarchy_keys(TRUE)[1], "srm_NrNotNAs") |>
  dplyr::distinct() |> dplyr::summarize(sum(srm_NrNotNAs)) |> dplyr::pull()
stopifnot(sum(!is.na(res[[config$get_response()[1]]])) == x)
res |> dplyr::select(c(config$hierarchy_keys(),"srm_NrNotNAs"  ,"srm_NrNotNARank")) |>
 dplyr::distinct() |>
 dplyr::arrange(!!!rlang::syms(c(config$hierarchy_keys()[1],"srm_NrNotNARank")))
#> # A tibble: 823 × 6
#>    protein_Id peptide_Id modPeptide_Id precursor_Id srm_NrNotNAs srm_NrNotNARank
#>    <chr>      <chr>      <chr>         <chr>               <int>           <int>
#>  1 A0A075B6I0 FSGSILGNK  _FSGSILGNK_   _FSGSILGNK_…           45               1
#>  2 A0A075B6P… DIVMTQSPL… _DIVMTQSPLSL… _DIVMTQSPLS…           45               1
#>  3 A0A075B6P… DIVMTQSPL… _DIVM[+16]TQ… _DIVM[+16]T…           45               1
#>  4 A0A075B6P… PGQSPQLLI… _PGQSPQLLIYL… _PGQSPQLLIY…           45               1
#>  5 A0A075B6P… PGQSPQLLI… _PGQSPQLLIYL… _PGQSPQLLIY…           45               1
#>  6 A0A075B6P… PGQSPQLLI… _PGQSPQLLIYL… _PGQSPQLLIY…           45               1
#>  7 A0A075B6P… SSQSLLHSN… _SSQSLLHSNGY… _SSQSLLHSNG…           45               1
#>  8 A0A075B6P… SSQSLLHSN… _SSQSLLHSNGY… _SSQSLLHSNG…           45               1
#>  9 A0A075B6R2 GLEWIGEIY… _GLEWIGEIYHS… _GLEWIGEIYH…           45               1
#> 10 A0A075B6R2 GLEWIGEIY… _GLEWIGEIYHS… _GLEWIGEIYH…           45               1
#> # ℹ 813 more rows
```
