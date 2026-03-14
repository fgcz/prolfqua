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
  config <- create_config_Skyline(isotopeLabel = "Isotope.Label.Type", ident_qValue = "Detection.Q.Value")
  config$factors[["Time"]] = "Sampling.Time.Point"
  return(config)
}
usethis::use_data(data_skylinePRMSample_A, overwrite = TRUE)

# --- data_skylineSRM_HL_A: fix config_f closure ---
load("data/data_skylineSRM_HL_A.rda")
data_skylineSRM_HL_A$config_f <- function() {
  skylineconfig_HL <- create_config_Skyline(isotopeLabel = "Isotope.Label", ident_qValue = "annotation_QValue")
  skylineconfig_HL$factors[["treatment_c"]] <- "Condition2"
  skylineconfig_HL$factors[["time_c"]] <- "time"
  skylineconfig_HL$is_response_transformed = FALSE
  return(skylineconfig_HL)
}
usethis::use_data(data_skylineSRM_HL_A, overwrite = TRUE)

# --- data_spectronautDIA250_A: fix config_f closure ---
load("data/data_spectronautDIA250_A.rda")
data_spectronautDIA250_A$config_f <- function() {
  spectronautDIAData250_config <- prolfqua::create_config_Spectronaut_Peptide(
    isotopeLabel = "Isotope.Label",
    ident_qValue = "EG.Qvalue"
  )
  spectronautDIAData250_config$factors[["coding"]] = "coding"
  spectronautDIAData250_config$factors[["sex"]] = "sex"
  spectronautDIAData250_config$factors[["age"]] = "age"
  spectronautDIAData250_config$factors[["Sample_id"]] = "Sample.Name"
  return(spectronautDIAData250_config)
}
usethis::use_data(data_spectronautDIA250_A, overwrite = TRUE)

cat("\nDone. All 5 data objects regenerated.\n")
