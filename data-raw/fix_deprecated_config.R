# Fix deprecated $table/$parameter active bindings in serialized data objects
#
# Root cause: During `R CMD INSTALL`, the lazyload DB creation phase evaluates
# all active bindings on R6 objects. The `$table` and `$parameter` active
# bindings on AnalysisConfiguration print deprecation messages when read.
#
# Sources of warnings:
#   1. data_ionstar.rda — IonstarData R6 with $config and $config_N (2 configs)
#   2. data_IonstarProtein_subsetNorm.rda — list with $config (1 config)
#   Total: 3 configs × 2 bindings = 6 deprecation messages
#
# Fix: Replace the noisy active bindings with silent versions that still return
# self for backwards compatibility. Uses unlockBinding() to allow replacement
# on R6's locked environments.
#
# Additionally fixes 3 .rda files with config_f closures that used the old
# config$table$factors[...] syntax.

library(prolfqua)

#' Replace deprecated active bindings on an AnalysisConfiguration with silent versions.
#' The R6 environment is locked, but unlockBinding() allows replacing individual bindings.
#' The new binding function uses the R6 enclosing environment so `self` resolves correctly.
silence_deprecated_bindings <- function(config) {
  stopifnot(inherits(config, "AnalysisConfiguration"))
  env <- config # R6 public env IS the object
  enclos <- env[[".__enclos_env__"]]
  for (nm in c("table", "parameter")) {
    if (exists(nm, envir = env, inherits = FALSE) && bindingIsActive(nm, env)) {
      unlockBinding(nm, env)
      silent_fn <- function() self
      environment(silent_fn) <- enclos
      makeActiveBinding(nm, silent_fn, env)
    }
  }
  invisible(config)
}

# --- data_ionstar: IonstarData R6 with config + config_N ---
load("data/data_ionstar.rda")
silence_deprecated_bindings(data_ionstar$config)
silence_deprecated_bindings(data_ionstar$config_N)
usethis::use_data(data_ionstar, overwrite = TRUE)

# --- data_IonstarProtein_subsetNorm: list with config ---
load("data/data_IonstarProtein_subsetNorm.rda")
silence_deprecated_bindings(data_IonstarProtein_subsetNorm$config)
usethis::use_data(data_IonstarProtein_subsetNorm, overwrite = TRUE)

# --- data_skylinePRMSample_A: fix config_f closure ---
load("data/data_skylinePRMSample_A.rda")
data_skylinePRMSample_A$config_f <- function() {
  config <- prolfqua::AnalysisConfiguration$new()
  config$fileName <- "Replicate.Name"
  config$hierarchy[["protein_Id"]] <- "Protein.Name"
  config$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  config$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
  config$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence", "Precursor.Charge", "Fragment.Ion", "Product.Charge")
  config$ident_qValue <- "Detection.Q.Value"
  config$set_response("Area")
  config$isotopeLabel <- "Isotope.Label.Type"
  config$factors[["Time"]] = "Sampling.Time.Point"
  return(config)
}
usethis::use_data(data_skylinePRMSample_A, overwrite = TRUE)

# --- data_skylineSRM_HL_A: fix config_f closure ---
load("data/data_skylineSRM_HL_A.rda")
data_skylineSRM_HL_A$config_f <- function() {
  config <- prolfqua::AnalysisConfiguration$new()
  config$fileName <- "Replicate.Name"
  config$hierarchy[["protein_Id"]] <- "Protein.Name"
  config$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  config$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
  config$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence", "Precursor.Charge", "Fragment.Ion", "Product.Charge")
  config$ident_qValue <- "annotation_QValue"
  config$set_response("Area")
  config$isotopeLabel <- "Isotope.Label"
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
  config$fileName <- "R.FileName"
  config$hierarchy[["protein_Id"]] <- "PG.ProteinAccessions"
  config$hierarchy[["peptide_Id"]] <- "PEP.StrippedSequence"
  config$hierarchy[["modPeptide_Id"]] <- "EG.ModifiedSequence"
  config$hierarchy[["precursor_Id"]] <- c("EG.ModifiedSequence", "FG.Charge")
  config$ident_qValue <- "EG.Qvalue"
  config$workIntensity <- "FG.Quantity"
  config$isotopeLabel <- "Isotope.Label"
  config$factors[["coding"]] = "coding"
  config$factors[["sex"]] = "sex"
  config$factors[["age"]] = "age"
  config$factors[["Sample_id"]] = "Sample.Name"
  return(config)
}
usethis::use_data(data_spectronautDIA250_A, overwrite = TRUE)

cat("\nDone. All 5 data objects regenerated.\n")
