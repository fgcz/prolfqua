context("test-tidyconfig_functions")

test_that("check config", {
  config <- create_config_Skyline()
  config$factors[["Time"]] = "Sampling.Time.Point"
  expect_equal(config$factor_keys(), "Time")
  expect_equal(config$hierarchy_keys(), c("protein_Id", "peptide_Id", "precursor_Id", "fragment_Id"))
})
