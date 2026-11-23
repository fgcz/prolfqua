# Reconstruct a configuration from minimal YAML

The inverse of
[`r6_extract_values()`](https://fgcz.github.io/prolfqua/reference/r6_extract_values.md):
rebuilds an
[AnalysisConfiguration](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md)
from the list read from a YAML file, in the current flat format or the
older nested one. Unknown keys raise an error.

## Usage

``` r
list_to_analysis_configuration(dd)
```

## Arguments

- dd:

  list with table and parameter elements as produced by
  r6_extract_values

## Value

an
[AnalysisConfiguration](https://fgcz.github.io/prolfqua/reference/AnalysisConfiguration.md)
with the fields set from `dd`.

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
configList <- prolfqua::r6_extract_values(DEAconfig)
stopifnot(class(configList) == "list")
config <- list_to_analysis_configuration(configList)
all.equal(prolfqua::r6_extract_values(config), configList)
#> [1] TRUE
```
