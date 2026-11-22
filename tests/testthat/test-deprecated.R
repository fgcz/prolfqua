## Renamed functions and simulator arguments (prolfqua 1.9.5) keep working with a deprecation warning

test_that("deprecated function names warn and return the same value as the new names", {
  sim <- sim_lfq_data_peptide_config(n_proteins = 5)
  lfq <- LFQData$new(sim$data, sim$config)
  config_list <- r6_extract_values(sim$config)

  expect_warning(old <- R6_extract_values(sim$config), "deprecated")
  expect_identical(old, config_list)
  expect_warning(old <- list_to_AnalysisConfiguration(config_list), "deprecated")
  expect_identical(r6_extract_values(old), r6_extract_values(list_to_analysis_configuration(config_list)))

  contrast_config <- ContrastConfiguration$new()
  contrast_list <- r6_extract_values(contrast_config)
  expect_warning(old <- list_to_ContrastConfiguration(contrast_list), "deprecated")
  expect_identical(r6_extract_values(old), r6_extract_values(list_to_contrast_configuration(contrast_list)))

  quantiles <- summarize_stats_quantiles(summarize_stats(lfq), lfq$relevant_factor_keys(), probs = c(0.5, 0.8))
  expect_warning(old <- lfq_power_t_test_quantiles_V2(quantiles$long), "deprecated")
  expect_identical(old, lfq_power_t_test_quantiles(quantiles$long))

  ropeca_input <- data.frame(
    contrast = "c1",
    protein_Id = rep(c("P1", "P2"), each = 3),
    estimate = c(1, 2, 3, -1, -2, 0),
    statistic = c(3, 1, 2, -0.5, -2, 0),
    p.value = c(0.01, 0.2, 0.03, 0.5, 0.04, 0.9),
    avgAbd = 1:6
  )
  args <- list(ropeca_input, subject_id = "protein_Id", estimate = "estimate", p.value = "p.value")
  expect_warning(old <- do.call(summary_ROPECA_median_p_scaled, args), "deprecated")
  expect_identical(old, do.call(summary_ropeca_median_p_scaled, args))

  expect_warning(p <- missigness_histogram(lfq), "deprecated")
  expect_s3_class(p, "ggplot")
  expect_warning(p <- lfq$get_Plotter()$missigness_histogram(), "deprecated")
  expect_s3_class(p, "ggplot")

  two_factor <- sim_lfq_data_2factor_config(n_proteins = 10, with_missing = FALSE)
  lfq2 <- LFQData$new(two_factor$data, two_factor$config)
  full <- build_model(lfq2, strategy_lm("abundance ~ Treatment + Background"))
  reduced <- build_model(lfq2, strategy_lm("abundance ~ Treatment"))
  args <- list(full$model_df, "full", reduced$model_df, "reduced")
  expect_warning(old <- do.call(LR_test, args), "lr_test")
  expect_identical(old, do.call(lr_test, args))

  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  expect_warning(old <- LFQDataToSummarizedExperiment(lfq), "lfqdata_to_summarized_experiment")
  expect_identical(old, lfqdata_to_summarized_experiment(lfq))
})

test_that("deprecated simulator arguments warn and give the same data as the new ones", {
  warnings <- testthat::capture_warnings(
    old <- withr::with_seed(1, sim_lfq_data(Nprot = 5, N = 3, PEPTIDE = TRUE))
  )
  expect_length(warnings, 3)
  expect_match(warnings, "deprecated")
  expect_identical(old, withr::with_seed(1, sim_lfq_data(n_proteins = 5, n_per_group = 3, peptide = TRUE)))

  expect_length(testthat::capture_warnings(old <- sim_lfq_data_peptide_config(Nprot = 5, N = 3)), 2)
  expect_identical(old$data, sim_lfq_data_peptide_config(n_proteins = 5, n_per_group = 3)$data)

  expect_warning(old <- sim_lfq_data_protein_config(Nprot = 5), "`Nprot` is deprecated")
  expect_identical(old$data, sim_lfq_data_protein_config(n_proteins = 5)$data)

  expect_length(
    testthat::capture_warnings(old <- sim_lfq_data_2factor_config(Nprot = 5, PEPTIDE = TRUE, TWO = FALSE)),
    3
  )
  expect_identical(old$data, sim_lfq_data_2factor_config(n_proteins = 5, peptide = TRUE, two_factor = FALSE)$data)
})

test_that("volcano functions take `data`; `.data` warns and gives the same plot", {
  d <- withr::with_seed(
    1,
    data.frame(
      fc = rnorm(40),
      BFDR = runif(40),
      condition = rep(c("A", "B"), 20),
      Prey = paste0("P", 1:40),
      modelName = "m",
      colour = "x",
      name = paste0("n", 1:40)
    )
  )
  expect_warning(old <- multigroup_volcano(.data = d, significance = "BFDR", label = "name"), "deprecated")
  new <- multigroup_volcano(data = d, significance = "BFDR", label = "name")
  expect_identical(ggplot2::layer_data(old, 1), ggplot2::layer_data(new, 1))

  skip_if_not_installed("plotly")
  skip_if_not_installed("crosstalk")
  expect_warning(old <- volcano_plotly(.data = d), "deprecated")
  new <- volcano_plotly(d)
  expect_length(old, length(new))
  expect_equal(plotly::plotly_build(old[[1]])$x$data[[1]]$x, plotly::plotly_build(new[[1]])$x$data[[1]]$x)
})

test_that("merge_contrasts_results warns on `model_name` and ignores it", {
  tab <- function(model, diff) {
    df <- data.frame(
      protein_Id = c("P1", "P2"),
      contrast = "A_vs_B",
      modelName = model,
      diff = diff,
      p.value = c(0.01, 0.2),
      FDR = c(0.02, 0.2),
      statistic = c(3, NA)
    )
    ContrastsTable$new(df, subject_id = "protein_Id", model_name = model)
  }
  prefer <- tab("prefer", c(1, 2))
  add <- tab("add", c(1.1, 2.1))
  expect_warning(old <- merge_contrasts_results(prefer, add, model_name = "my_merge"), "deprecated")
  new <- merge_contrasts_results(prefer, add)
  expect_identical(old$merged$get_contrasts(), new$merged$get_contrasts())
  expect_identical(old$merged$model_name, "prefer_add")
})
