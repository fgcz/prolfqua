# Rebuild the configurations of data_ionstar and add identification confidence columns
#
# Problem: data_ionstar holds configurations serialized by old prolfqua versions (camelCase fields,
# `ident_q_value = "pep"`). prolfqua renamed `ident_q_value` (smaller is better) to `identification_confidence`
# (in [0, 1], larger is better, used as model weight) and added the feature-level columns
# `identification_confidence_feature` and `nr_children_feature`.
#
# Fix: build fresh configurations from the old field values, set identification_confidence = 1 - pep (pep is
# missing only where the intensity is missing), and derive the feature-level columns as the maximum over samples.
# Run from the package root after data-raw/data_ionstar_nr_children.R.

devtools::load_all(".")
load("data/data_ionstar.rda")

rebuild_config <- function(old) {
  config <- AnalysisConfiguration$new()
  config$file_name <- old$file_name
  config$sample_name <- old$sample_name
  config$hierarchy <- old$hierarchy
  config$hierarchy_depth <- old$hierarchy_depth
  config$factors <- old$factors
  config$factor_depth <- old$factor_depth
  config$work_intensity <- old$work_intensity
  config$is_response_transformed <- old$is_response_transformed
  config$nr_children <- old$nr_children
  config
}

max_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

add_confidence <- function(data, config) {
  data$identification_confidence <- 1 - data$pep
  dplyr::mutate(
    data,
    identification_confidence_feature = max_or_na(.data$identification_confidence),
    nr_children_feature = max_or_na(.data$nr_children),
    .by = dplyr::all_of(config$hierarchy_keys())
  )
}

data_ionstar$config <- rebuild_config(data_ionstar$config)
data_ionstar$config_N <- rebuild_config(data_ionstar$config_N)
data_ionstar$data <- add_confidence(data_ionstar$data, data_ionstar$config)
data_ionstar$data_N <- add_confidence(data_ionstar$data_N, data_ionstar$config_N)
save(data_ionstar, file = "data/data_ionstar.rda", compress = "xz")
