test_that("missingness_histogram facets by the factors it is given (D-review fold-in)", {
  sim <- sim_lfq_data_2factor_config(n_proteins = 10)
  lfq <- LFQData$new(sim$data, sim$config)
  p <- .missingness_histogram(lfq, factors = "Treatment")
  facet_vars <- c(names(p$facet$params$rows), names(p$facet$params$cols))
  expect_true("Treatment" %in% facet_vars)
  expect_false("Background" %in% facet_vars)
  expect_true(all(facet_vars %in% colnames(p$data)))
})
