# Convert old proflqua configurations (prolfqua 0.4) to new Analysis configurations prolfqua 0.5

Convert old proflqua configurations (prolfqua 0.4) to new Analysis
configurations prolfqua 0.5

## Usage

``` r
old2new(config)
```

## Arguments

- config:

  an old (as stored in the example data) AnalysisConfiguration

## Value

a new AnalysisConfiguration

## Examples

``` r
dd <- prolfqua_data("data_ionstar")$filtered()
#> Column added : nr_peptide_Id_IN_protein_Id
dd$config <- old2new(dd$config)
```
