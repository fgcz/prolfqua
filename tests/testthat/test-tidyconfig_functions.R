context("test-tidyconfig_functions")

test_that("check config", {
  config <- createSkylineConfiguration()
  config$table$factors[["Time"]] = "Sampling.Time.Point"
  expect_equal(config$table$factorKeys(),"Time")
  expect_equal(config$table$hierarchyKeys(),c("protein_Id","peptide_Id","precursor_Id","fragment_Id"))
})
