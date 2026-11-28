# Rebuild the configurations of data_ionstar with the peptide consistency fields
#
# Problem: data_ionstar holds AnalysisConfiguration objects serialized before the fields `weight_columns`,
# `consistency_factor` and `use_peptide_consistency` existed. A serialized R6 object keeps its old set of fields and
# is locked, so the aggregations, which set `weight_columns` and `consistency_factor`, fail on it.
#
# Fix: rebuild each configuration from its values with the current class definition. The data are unchanged.
# Run from the package root.

devtools::load_all(".")
load("data/data_ionstar.rda")

rebuild_config <- function(old) {
  list_to_analysis_configuration(r6_extract_values(old))
}

data_ionstar$config <- rebuild_config(data_ionstar$config)
data_ionstar$config_N <- rebuild_config(data_ionstar$config_N)
save(data_ionstar, file = "data/data_ionstar.rda", compress = "xz")
