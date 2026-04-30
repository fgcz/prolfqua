# Setup a tidy table compatible with a [`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md)

Extracts columns relevant for a configuration from a data frame and
create new columns e.g. sampleName column etc.

## Usage

``` r
setup_analysis(data, configuration, cc = TRUE, from_factors = FALSE)
```

## Arguments

- data:

  data.frame

- configuration:

  AnalysisConfiguration

- cc:

  complete cases default TRUE

- from_factors:

  if TRUE, create sampleName from factor columns

## Value

The computed result.

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
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Examples

``` r
skylineconfig <- AnalysisConfiguration$new()
skylineconfig$file_name <- "Replicate.Name"
skylineconfig$hierarchy[["protein_Id"]] <- "Protein.Name"
skylineconfig$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
skylineconfig$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
skylineconfig$hierarchy[["fragment_Id"]] <- c(
  "Peptide.Sequence", "Precursor.Charge",
  "Fragment.Ion", "Product.Charge"
)
skylineconfig$ident_q_value <- "Detection.Q.Value"
skylineconfig$set_response("Area")
skylineconfig$isotope_label <- "Isotope.Label.Type"
skylineconfig$factors[["Time"]] = "Sampling.Time.Point"
sample_analysis <- setup_analysis(prolfqua_data('data_skylinePRMSample_A')$data, skylineconfig)
#> creating sampleName from file_name column
#> Warning: no nr_children column specified in the data, adding column nr_children and setting to 1.
#> completing cases
#> completing cases done
#> setup done

# Example with normValue column (e.g., creatinine)
set.seed(1234)
data <- sim_lfq_data(Nprot = 10, PEPTIDE = TRUE, N = 4)
data$nr_children <- 1
data$isotopeLabel <- "light"
data$qValue <- 0

# Add creatinine values per sample (not per protein/peptide)
sample_creatinine <- data |> dplyr::select(sample) |> dplyr::distinct() |>
  dplyr::mutate(creatinine = rnorm(dplyr::n(), mean = 100, sd = 10))
data <- dplyr::inner_join(data, sample_creatinine, by = "sample")

config <- AnalysisConfiguration$new()
config$file_name = "sample"
config$factors["group_"] = "group"
config$hierarchy[["protein_Id"]] = c("proteinID", "idtype2")
config$hierarchy[["peptide_Id"]] = "peptideID"
config$set_response("abundance")
config$norm_value = "creatinine"

adata <- setup_analysis(data, config)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
```
