# Generate instances of AnalysisConfiguration

configurations examples of or various signal processing software outputs

file must be read with tidyMQ_Peptides, you will still need to add the
factors (explanatory variables).

## Usage

``` r
create_config_Skyline(
  isotopeLabel = "Isotope.Label",
  ident_qValue = "annotation_QValue"
)

create_config_Spectronaut_Peptide(
  isotopeLabel = "Isotope.Label",
  ident_qValue = "EG.Qvalue"
)

create_config_MQ_peptide(
  ident_qValue = "pep",
  intensity = "peptide.intensity",
  isotopeLabel = "isotope"
)

create_config_MSFragger_MSstats()
```

## Arguments

- isotopeLabel:

  isotope

- ident_qValue:

  pep

- intensity:

  peptide.intensity

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_factors()`](https://wolski.github.io/prolfqua/reference/separate_factors.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`spread_response_by_IsotopeLabel()`](https://wolski.github.io/prolfqua/reference/spread_response_by_IsotopeLabel.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
skylineconfig <- create_config_Skyline()
skylineconfig$factors[["Time"]] = "Sampling.Time.Point"
skylineconfig$factor_keys()
#> [1] "Time"
skylineconfig$hierarchy_keys()
#> [1] "protein_Id"   "peptide_Id"   "precursor_Id" "fragment_Id" 



spectronautconfig <- create_config_Spectronaut_Peptide()
config <- create_config_Spectronaut_Peptide()
config$factors[["coding"]] = "coding"
config$factors[["sex"]] = "sex"
config$factors[["age"]] = "age"
config$factors[["Sample_id"]] = "Sample.Name"

tmp <- create_config_MQ_peptide()

create_config_MSFragger_MSstats()
#> <AnalysisConfiguration>
#>   Public:
#>     annotation_vars: function () 
#>     bin_resp: 
#>     clone: function (deep = FALSE) 
#>     factorDepth: 1
#>     factor_keys: function () 
#>     factor_keys_depth: function () 
#>     factors: list
#>     fileName: Run
#>     get_response: function () 
#>     hierarchy: list
#>     hierarchyDepth: 1
#>     hierarchyKeys: function (rev = FALSE) 
#>     hierarchy_keys: function (rev = FALSE) 
#>     hierarchy_keys_depth: function (names = TRUE) 
#>     hkeysDepth: function (names = TRUE) 
#>     id_required: function () 
#>     id_vars: function () 
#>     ident_Score: 
#>     ident_qValue: pep
#>     initialize: function (analysisTableAnnotation = NULL, analysisParameter = NULL) 
#>     is_response_transformed: FALSE
#>     isotopeLabel: IsotopeLabelType
#>     min_peptides_protein: 2
#>     normValue: NULL
#>     nr_children: nr_children
#>     opt_mz: 
#>     opt_rt: 
#>     parameter: active binding
#>     pop_response: function () 
#>     sampleName: sampleName
#>     sep: ~
#>     set_response: function (colName) 
#>     table: active binding
#>     value_vars: function () 
#>     workIntensity: Intensity
```
