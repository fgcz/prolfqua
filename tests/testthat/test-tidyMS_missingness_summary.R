test_that("missingness_histogram facets by the factors it is given (D-review fold-in)", {
  sim <- sim_lfq_data_2factor_config(n_proteins = 10)
  lfq <- LFQData$new(sim$data, sim$config)
  p <- .missingness_histogram(lfq, factors = "Treatment")
  facet_vars <- c(names(p$facet$params$rows), names(p$facet$params$cols))
  expect_true("Treatment" %in% facet_vars)
  expect_false("Background" %in% facet_vars)
  expect_true(all(facet_vars %in% colnames(p$data)))
})

test_that("summariser missingness wrappers return their underlying tables", {
  sim <- sim_lfq_data_protein_config(n_proteins = 10, weight_missing = 0.4, seed = 1)
  lfq <- LFQData$new(sim$data, sim$config)
  summariser <- lfq$get_Summariser()
  expect_identical(summariser$interaction_missing_stats(), summarize_stats_factors(lfq))
  expect_identical(summariser$missingness_per_group(), .missingness_per_condition(lfq)$data)
  expect_identical(summariser$missingness_per_group_cumsum(), .missingness_per_condition_cumsum(lfq)$data)
})
