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
[`PACKAGE_DATA`](https://fgcz.github.io/prolfqua/reference/PACKAGE_DATA.md),
[`data_checksummarizationrobust87`](https://fgcz.github.io/prolfqua/reference/data_checksummarizationrobust87.md),
[`data_checksummarizerobust69`](https://fgcz.github.io/prolfqua/reference/data_checksummarizerobust69.md),
[`data_ionstar`](https://fgcz.github.io/prolfqua/reference/data_ionstar.md),
[`x5463yzwer453bbb`](https://fgcz.github.io/prolfqua/reference/x5463yzwer453bbb.md)

## Examples

``` r
ionstar <- prolfqua_data("data_ionstar")
names(ionstar)
#>  [1] ".__enclos_env__"   "subset_normalized" "normalized"       
#>  [4] "data_N"            "clone"             "initialize"       
#>  [7] "config"            "data"              "Pep"              
#> [10] "config_N"          "filtered"         
```
