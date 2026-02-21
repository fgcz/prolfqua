# Extract response column of a protein into matrix

Used to apply the median polish function working on matrices to a tidy
table

## Usage

``` r
response_as_matrix(pdata, config)
```

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following object is masked from ‘package:MASS’:
#> 
#>     select
#> The following object is masked from ‘package:Biobase’:
#> 
#>     combine
#> The following objects are masked from ‘package:GenomicRanges’:
#> 
#>     intersect, setdiff, union
#> The following object is masked from ‘package:Seqinfo’:
#> 
#>     intersect
#> The following objects are masked from ‘package:IRanges’:
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from ‘package:S4Vectors’:
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from ‘package:BiocGenerics’:
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from ‘package:generics’:
#> 
#>     explain
#> The following object is masked from ‘package:matrixStats’:
#> 
#>     count
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

bb <- prolfqua_data("data_ionstar")$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
bb$config <- old2new(bb$config)
stopifnot(nrow(bb$data) == 25780)
configur <- bb$config
data <- bb$data

xnested <- data |>
  dplyr::group_by(across(all_of(configur$hierarchy_keys_depth()))) |>
  tidyr::nest()
x <- xnested$data[[1]]
nn <- x |>
  dplyr::select(base::setdiff(
    configur$hierarchy_keys(),
    configur$hierarchy_keys_depth()
  )) |>
  dplyr::distinct() |>
  nrow()

xx <- response_as_matrix(x, configur)
stopifnot(dim(xx) == c(nn, 20))

# change hierarchyDepth ###################
conf <- configur$clone(deep = TRUE)
conf$hierarchyDepth <- 1

xnested <- data |>
  dplyr::group_by(across(all_of(conf$hierarchy_keys_depth()))) |>
  tidyr::nest()

x <- xnested$data[[1]]
nn <- x |>
  dplyr::select(base::setdiff(
    configur$hierarchy_keys(),
    configur$hierarchy_keys_depth()
  )) |>
  dplyr::distinct() |>
  nrow()

xx <- response_as_matrix(x, conf)
stopifnot(dim(xx) == c(nn, 20))
```
