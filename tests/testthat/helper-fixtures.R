# Shared test fixtures.

# Protein-level LFQData from sim_lfq_data_protein_config(), response renamed to transformedIntensity;
# `...` goes to the simulator (e.g. weight_missing).
protein_lfqdata_fixture <- function(n_proteins = 30, seed = 42, ...) {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = n_proteins, seed = seed, ...)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  lfqdata
}

# Analysis of a 20-protein peptide simulation with `method`, given the input shape the method needs
facade_method_case <- function(method) {
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, seed = 42)
  peptide <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  if (method == "limpa") {
    input <- AggregateLimpa$new(peptide, "protein")$aggregate()
  } else if (lookup_facade(method)$needs == "nested") {
    input <- peptide
  } else {
    input <- peptide$get_Aggregator()$aggregate()
  }
  model <- if (method == "lmer_nested") "~ group_ + (1 | peptide_Id) + (1 | sampleName)" else "~ group_"
  build_contrast_analysis(input, model, c(A_vs_Ctrl = "group_A - group_Ctrl"), method)
}

# Suggests packages each method of facade_method_case() needs
facade_method_dependencies <- list(
  rfit = "Rfit",
  rfit_impute = "Rfit",
  lmer_nested = c("lme4", "lmerTest"),
  firth = "logistf",
  limpa = "limpa",
  limpa_nested = "limpa"
)
