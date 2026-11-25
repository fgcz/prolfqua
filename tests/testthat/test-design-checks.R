# Design checks before fitting: collinear designs are rejected for every method (srr RE2.4, RE7.0).

#' @srrstats {RE7.0, RE7.0a} Tests build designs with exact relationships between predictors (a factor duplicating
#'   another, an interaction cell without samples) and confirm that every method rejects them with
#'   `prolfqua_error_collinear_design` before fitting.
#' @noRd
NULL

aliased_data <- function(simulate) {
  sim <- simulate(n_proteins = 5, with_missing = FALSE, seed = 42)
  data <- sim$data
  data$batch <- c(A = "x", B = "y", Ctrl = "z")[data$group_]
  sim$config$factors[["batch"]] <- "batch"
  suppressMessages(LFQData$new(data, sim$config)$get_Transformer()$log2()$lfq)
}

test_that("a factor aliased with another is rejected by every built-in method", {
  withr::local_options(prolfqua.verbose = FALSE)
  data <- list(
    aggregated = aliased_data(sim_lfq_data_protein_config),
    nested = aliased_data(sim_lfq_data_peptide_config)
  )
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  builtin <- Filter(function(facade) facade$package == "prolfqua", list_facades())
  # the design check runs before any backend package is needed
  for (method in names(builtin)) {
    expect_error(
      build_contrast_analysis(data[[builtin[[method]]$needs]], "~ group_ + batch", contrasts, method),
      "`batchy`, `batchz` are linear combinations",
      class = "prolfqua_error_collinear_design",
      label = method
    )
  }
  expect_no_error(suppressWarnings(build_contrast_analysis(data$aggregated, "~ group_", contrasts, "lm")))
})

test_that("nested methods check the fixed effects only", {
  withr::local_options(prolfqua.verbose = FALSE)
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  peptide <- aliased_data(sim_lfq_data_peptide_config)
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  expect_error(
    build_contrast_analysis(peptide, "~ group_ + batch + (1 | peptide_Id)", contrasts, "lmer_nested"),
    class = "prolfqua_error_collinear_design"
  )
  fit <- suppressMessages(suppressWarnings(
    build_contrast_analysis(peptide, "~ group_ + (1 | peptide_Id)", contrasts, "lmer_nested")
  ))
  expect_s3_class(fit$get_contrasts(), "data.frame")
})

test_that("an interaction cell without samples is rejected; a complete two-factor design is accepted", {
  withr::local_options(prolfqua.verbose = FALSE)
  sim <- sim_lfq_data_2factor_config(n_proteins = 5, with_missing = FALSE)
  lfq <- LFQData$new(sim$data, sim$config)
  contrasts <- c(TB_vs_A = "TreatmentB - TreatmentA")
  expect_no_error(build_contrast_analysis(lfq, "~ Treatment * Background", contrasts, "lm"))
  empty_cell <- lfq$get_copy()
  empty_cell$set_data(dplyr::filter(lfq$data_long(), !(Treatment == "B" & Background == "Z")))
  expect_error(
    build_contrast_analysis(empty_cell, "~ Treatment * Background", contrasts, "lm"),
    "TreatmentB:BackgroundZ",
    class = "prolfqua_error_collinear_design"
  )
  expect_no_error(build_contrast_analysis(empty_cell, "~ Treatment + Background", contrasts, "lm"))
})
