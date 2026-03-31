# Re-serialize data objects with current AnalysisConfiguration class definition
#
# Problem: data_ionstar and data_IonstarProtein_subsetNorm contain serialized
# R6 objects with the OLD AnalysisConfiguration class definition (camelCase
# field names). When loaded, these stale R6 objects lack the new snake_case
# fields and active binding aliases.
#
# Fix: Extract config values via R6_extract_values(), rebuild with
# list_to_AnalysisConfiguration(), and re-save.
#
# The config_f closure datasets (skyline*, spectronaut*) are fine — they
# create fresh AnalysisConfiguration objects at call time.

library(prolfqua)

#' Rebuild a stale serialized AnalysisConfiguration using current class definition
rebuild_config <- function(config) {
  stopifnot(inherits(config, "AnalysisConfiguration"))
  vals <- R6_extract_values(config)
  list_to_AnalysisConfiguration(vals)
}

# --- data_ionstar: IonstarData R6 with config + config_N ---
load("data/data_ionstar.rda")
data_ionstar$config <- rebuild_config(data_ionstar$config)
data_ionstar$config_N <- rebuild_config(data_ionstar$config_N)
usethis::use_data(data_ionstar, overwrite = TRUE)

# --- data_IonstarProtein_subsetNorm: list with config ---
load("data/data_IonstarProtein_subsetNorm.rda")
data_IonstarProtein_subsetNorm$config <- rebuild_config(data_IonstarProtein_subsetNorm$config)
usethis::use_data(data_IonstarProtein_subsetNorm, overwrite = TRUE)

# --- data_skylinePRMSample_A: fix config_f closure ---
load("data/data_skylinePRMSample_A.rda")
data_skylinePRMSample_A$config_f <- function() {
  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "Replicate.Name"
  config$hierarchy[["protein_Id"]] <- "Protein.Name"
  config$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  config$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
  config$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence", "Precursor.Charge", "Fragment.Ion", "Product.Charge")
  config$ident_q_value <- "Detection.Q.Value"
  config$set_response("Area")
  config$isotope_label <- "Isotope.Label.Type"
  config$factors[["Time"]] = "Sampling.Time.Point"
  return(config)
}
usethis::use_data(data_skylinePRMSample_A, overwrite = TRUE)

# --- data_skylineSRM_HL_A: fix config_f closure ---
load("data/data_skylineSRM_HL_A.rda")
data_skylineSRM_HL_A$config_f <- function() {
  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "Replicate.Name"
  config$hierarchy[["protein_Id"]] <- "Protein.Name"
  config$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  config$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
  config$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence", "Precursor.Charge", "Fragment.Ion", "Product.Charge")
  config$ident_q_value <- "annotation_QValue"
  config$set_response("Area")
  config$isotope_label <- "Isotope.Label"
  config$factors[["treatment_c"]] <- "Condition2"
  config$factors[["time_c"]] <- "time"
  config$is_response_transformed = FALSE
  return(config)
}
usethis::use_data(data_skylineSRM_HL_A, overwrite = TRUE)

# --- data_spectronautDIA250_A: fix config_f closure ---
load("data/data_spectronautDIA250_A.rda")
data_spectronautDIA250_A$config_f <- function() {
  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "R.FileName"
  config$hierarchy[["protein_Id"]] <- "PG.ProteinAccessions"
  config$hierarchy[["peptide_Id"]] <- "PEP.StrippedSequence"
  config$hierarchy[["modPeptide_Id"]] <- "EG.ModifiedSequence"
  config$hierarchy[["precursor_Id"]] <- c("EG.ModifiedSequence", "FG.Charge")
  config$ident_q_value <- "EG.Qvalue"
  config$work_intensity <- "FG.Quantity"
  config$isotope_label <- "Isotope.Label"
  config$factors[["coding"]] = "coding"
  config$factors[["sex"]] = "sex"
  config$factors[["age"]] = "age"
  config$factors[["Sample_id"]] = "Sample.Name"
  return(config)
}
usethis::use_data(data_spectronautDIA250_A, overwrite = TRUE)

cat("\nDone. All 5 data objects regenerated.\n")
