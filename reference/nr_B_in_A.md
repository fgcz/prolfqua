# Compute nr of B per A

Compute nr of B per A

## Usage

``` r
nr_B_in_A(pdata, config, merge = TRUE)
```

## Arguments

- pdata:

  data.frame

- config:

  AnalysisConfiguration

## Examples

``` r
bb <- sim_lfq_data_peptide_config(Nprot = 100)
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- bb$config$clone(deep=TRUE)
data <- bb$data
hierarchy <- config$hierarchy_keys()
res <- nr_B_in_A(data, config)
#> Column added : nr_peptide_Id_IN_protein_Id

res$data |>
  dplyr::select(all_of(c(config$hierarchy_keys_depth(),  res$name))) |>
  dplyr::distinct() |>
  dplyr::pull() |> table()
#> 
#>  1  2  3  4  5  6  7  8  9 10 12 
#> 22 26 16 10  8  2  7  4  1  1  3 


bb <- prolfqua::prolfqua_data('data_skylineSRM_HL_A')
config <- old2new(bb$config_f())
#> config$table is deprecated, use config directly
#> config$table is deprecated, use config directly
#> config$table is deprecated, use config directly
data <- bb$data
data$Area[data$Area == 0] <- NA
analysis <- setup_analysis(data, config)
#> creating sampleName from fileName column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done

resDataStart <- bb$analysis(bb$data, config)
#> creating sampleName from fileName column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done

nr_B_in_A(resDataStart, config)
#> Column added : nr_peptide_Id_IN_protein_Id
#> $data
#> # A tibble: 5,715 × 13
#>    Replicate.Name         sampleName treatment_c time_c Isotope.Label protein_Id
#>    <chr>                  <chr>      <chr>       <chr>  <chr>         <chr>     
#>  1 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#>  2 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#>  3 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#>  4 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#>  5 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#>  6 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#>  7 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#>  8 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#>  9 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#> 10 04_S172853_lowG_p06_d… 04_S17285… lowG        p06    heavy         sp|O95477…
#> # ℹ 5,705 more rows
#> # ℹ 7 more variables: peptide_Id <chr>, precursor_Id <chr>, fragment_Id <chr>,
#> #   Area <dbl>, annotation_QValue <dbl>, nr_children <dbl>,
#> #   nr_peptide_Id_IN_protein_Id <int>
#> 
#> $name
#> [1] "nr_peptide_Id_IN_protein_Id"
#> 
nr_B_in_A(resDataStart, config, merge = FALSE)
#> # A tibble: 13 × 2
#>    protein_Id            nr_peptide_Id_IN_protein_Id
#>    <chr>                                       <int>
#>  1 iRT_Protein                                    11
#>  2 sp|O95477|ABCA1_HUMAN                           3
#>  3 sp|P02786|TFR1_HUMAN                            3
#>  4 sp|P04406|G3P_HUMAN                             3
#>  5 sp|P08195|4F2_HUMAN                             3
#>  6 sp|P49281|NRAM2_HUMAN                           3
#>  7 sp|P63104|1433Z_HUMAN                           2
#>  8 sp|Q01650|LAT1_HUMAN                            3
#>  9 sp|Q15043|S39AE_HUMAN                           3
#> 10 sp|Q15365|PCBP1_HUMAN                           3
#> 11 sp|Q9C0K1|S39A8_HUMAN                           2
#> 12 sp|Q9NP59|S40A1_HUMAN                           3
#> 13 sp|Q9UHI5|LAT2_HUMAN                            3
config$hierarchyDepth <- 2
nr_B_in_A(resDataStart, config, merge = FALSE)
#> # A tibble: 45 × 3
#> # Groups:   protein_Id [13]
#>    protein_Id  peptide_Id     nr_precursor_Id_IN_protein_Id_peptide_Id
#>    <chr>       <chr>                                             <int>
#>  1 iRT_Protein ADVTPADFSEWSK                                         1
#>  2 iRT_Protein DGLDAASYYAPVR                                         1
#>  3 iRT_Protein GAGSSEPVTGLDAK                                        1
#>  4 iRT_Protein GTFIIDPAAVIR                                          1
#>  5 iRT_Protein GTFIIDPGGVIR                                          1
#>  6 iRT_Protein LFLQFGAQGSPFLK                                        1
#>  7 iRT_Protein LGGNEQVTR                                             1
#>  8 iRT_Protein TPVISGGPYEYR                                          1
#>  9 iRT_Protein TPVITGAPYEYR                                          1
#> 10 iRT_Protein VEATFGVDESNAK                                         1
#> # ℹ 35 more rows
```
