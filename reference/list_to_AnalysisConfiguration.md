# read minimal yaml to reconstruct configuration

read minimal yaml to reconstruct configuration

## Usage

``` r
list_to_AnalysisConfiguration(dd)
```

## Arguments

- dd:

  list with table and parameter elements as produced by
  R6_extract_values

## Examples

``` r
DEAconfig <- AnalysisConfiguration$new()
DEAconfig$fileName <- "Replicate.Name"
DEAconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
DEAconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
DEAconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
DEAconfig$hierarchy[["fragment_Id"]] <- c(
  "Peptide.Sequence", "Precursor.Charge",
  "Fragment.Ion", "Product.Charge"
)
DEAconfig$ident_qValue <- "annotation_QValue"
DEAconfig$set_response("Area")
DEAconfig$isotopeLabel <- "Isotope.Label"
configList <- prolfqua::R6_extract_values(DEAconfig)
#> config$parameter is deprecated, use config directly
#> config$table is deprecated, use config directly
#> config$parameter is deprecated, use config directly
#> config$table is deprecated, use config directly
stopifnot(class(configList) == "list")
config <- list_to_AnalysisConfiguration(configList)
all.equal(prolfqua::R6_extract_values(config), configList)
#> config$parameter is deprecated, use config directly
#> config$table is deprecated, use config directly
#> config$parameter is deprecated, use config directly
#> config$table is deprecated, use config directly
#> [1] TRUE
```
