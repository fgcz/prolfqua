# load data from prolfqua

load data from prolfqua

## Usage

``` r
prolfqua_data(datastr, package = "prolfqua")
```

## Arguments

- datastr:

  name of dataset

- package:

  default prolfqua

## Value

The computed result.

## See also

Other data:
[`PACKAGE_DATA`](https://wolski.github.io/prolfqua/reference/PACKAGE_DATA.md),
[`data_IonstarProtein_subsetNorm`](https://wolski.github.io/prolfqua/reference/data_IonstarProtein_subsetNorm.md),
[`data_SAINTe_output`](https://wolski.github.io/prolfqua/reference/data_SAINTe_output.md),
[`data_benchmarkExample`](https://wolski.github.io/prolfqua/reference/data_benchmarkExample.md),
[`data_checksummarizationrobust87`](https://wolski.github.io/prolfqua/reference/data_checksummarizationrobust87.md),
[`data_checksummarizerobust`](https://wolski.github.io/prolfqua/reference/data_checksummarizerobust.md),
[`data_checksummarizerobust69`](https://wolski.github.io/prolfqua/reference/data_checksummarizerobust69.md),
[`data_correlatedPeptideList`](https://wolski.github.io/prolfqua/reference/data_correlatedPeptideList.md),
[`data_ionstar`](https://wolski.github.io/prolfqua/reference/data_ionstar.md),
[`data_skylinePRMSample_A`](https://wolski.github.io/prolfqua/reference/data_skylinePRMSample_A.md),
[`data_skylineSRM_HL_A`](https://wolski.github.io/prolfqua/reference/data_skylineSRM_HL_A.md),
[`data_spectronautDIA250_A`](https://wolski.github.io/prolfqua/reference/data_spectronautDIA250_A.md),
[`data_test_confusion_matrix_scores`](https://wolski.github.io/prolfqua/reference/data_test_confusion_matrix_scores.md),
[`x5463yzwer453bbb`](https://wolski.github.io/prolfqua/reference/x5463yzwer453bbb.md)

## Examples

``` r
ionstar <- prolfqua_data("data_ionstar")
names(ionstar)
#>  [1] ".__enclos_env__"   "subset_normalized" "normalized"       
#>  [4] "data_N"            "clone"             "initialize"       
#>  [7] "config"            "data"              "Pep"              
#> [10] "config_N"          "filtered"         
```
