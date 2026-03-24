# Generate instances of AnalysisConfiguration

Configuration examples for various signal processing software outputs.

file must be read with tidyMQ_Peptides, you will still need to add the
factors (explanatory variables).

## Usage

``` r
create_config_MQ_peptide(
  ident_qValue = "pep",
  intensity = "peptide.intensity",
  isotopeLabel = "isotope"
)
```

## Arguments

- ident_qValue:

  pep

- intensity:

  peptide.intensity

- isotopeLabel:

  isotope

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
# Skyline configuration
skylineconfig <- AnalysisConfiguration$new()
skylineconfig$fileName <- "Replicate.Name"
skylineconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
skylineconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
skylineconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
skylineconfig$hierarchy[["fragment_Id"]] <- c(
  "Peptide.Sequence", "Precursor.Charge",
  "Fragment.Ion", "Product.Charge"
)
skylineconfig$ident_qValue <- "annotation_QValue"
skylineconfig$set_response("Area")
skylineconfig$isotopeLabel <- "Isotope.Label"
skylineconfig$factors[["Time"]] = "Sampling.Time.Point"
skylineconfig$factor_keys()
#> [1] "Time"
skylineconfig$hierarchy_keys()
#> [1] "protein_Id"   "peptide_Id"   "precursor_Id" "fragment_Id" 

# Spectronaut configuration
spectronautconfig <- AnalysisConfiguration$new()
spectronautconfig$fileName <- "R.FileName"
spectronautconfig$hierarchy[["protein_Id"]] <- "PG.ProteinAccessions"
spectronautconfig$hierarchy[["peptide_Id"]] <- "PEP.StrippedSequence"
spectronautconfig$hierarchy[["modPeptide_Id"]] <- "EG.ModifiedSequence"
spectronautconfig$hierarchy[["precursor_Id"]] <- c("EG.ModifiedSequence", "FG.Charge")
spectronautconfig$ident_qValue <- "EG.Qvalue"
spectronautconfig$workIntensity <- "FG.Quantity"
spectronautconfig$isotopeLabel <- "Isotope.Label"
spectronautconfig$factors[["coding"]] = "coding"
spectronautconfig$factors[["sex"]] = "sex"
spectronautconfig$factors[["age"]] = "age"
spectronautconfig$factors[["Sample_id"]] = "Sample.Name"
tmp <- create_config_MQ_peptide()
```
