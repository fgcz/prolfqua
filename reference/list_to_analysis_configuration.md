# Reconstruct a configuration from minimal YAML

Reconstruct a configuration from minimal YAML

## Usage

``` r
list_to_analysis_configuration(dd)
```

## Arguments

- dd:

  list with table and parameter elements as produced by
  r6_extract_values

## Value

The computed result.

## Examples

``` r
DEAconfig <- AnalysisConfiguration$new()
DEAconfig$file_name <- "Replicate.Name"
DEAconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
DEAconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
DEAconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
DEAconfig$hierarchy[["fragment_Id"]] <- c(
  "Peptide.Sequence", "Precursor.Charge",
  "Fragment.Ion", "Product.Charge"
)
DEAconfig$ident_q_value <- "annotation_QValue"
DEAconfig$set_response("Area")
DEAconfig$isotope_label <- "Isotope.Label"
configList <- prolfqua::r6_extract_values(DEAconfig)
stopifnot(class(configList) == "list")
config <- list_to_analysis_configuration(configList)
all.equal(prolfqua::r6_extract_values(config), configList)
#> [1] TRUE
```
