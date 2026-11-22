# Shared test fixtures.

# Protein-level LFQData from sim_lfq_data_protein_config(), response renamed to transformedIntensity;
# `...` goes to the simulator (e.g. weight_missing).
protein_lfqdata_fixture <- function(n_proteins = 30, seed = 42, ...) {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = n_proteins, seed = seed, ...)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  lfqdata
}
