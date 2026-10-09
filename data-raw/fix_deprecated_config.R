# Re-serialize data objects with current AnalysisConfiguration class definition
#
# Problem: data_ionstar contains serialized R6 objects with the OLD
# AnalysisConfiguration class definition (camelCase field names). When loaded,
# these stale R6 objects lack the new snake_case fields and active binding
# aliases.
#
# Fix: Extract config values via R6_extract_values(), rebuild with
# list_to_AnalysisConfiguration(), and re-save.
#
# data_IonstarProtein_subsetNorm moved to prolfquasaint in prolfqua 1.9.0; the
# skyline* and spectronaut* datasets were removed.

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

cat("\nDone. data_ionstar regenerated.\n")
