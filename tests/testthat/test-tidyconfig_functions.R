context("test-tidyconfig_functions")

test_that("check config", {
  config <- createSkylineConfiguration()
  config$table$factors[["Time"]] = "Sampling.Time.Point"
  expect_equal(config$table$factorKeys(),"Time")
  expect_equal(config$table$hierarchyKeys(),c("protein_Id","peptide_Id","precursor_Id","fragment_Id"))
})

test_that("run test_resultsV1_2954_modelling.R",{
  source("script_resultsV1_2954_modelling.R", echo=FALSE)
})

test_that("run script_resultsV1_2954_create.R",{
  source("tests/testthat/script_resultsV1_2954_create.R", echo=FALSE)
})
