context("test-tidyconfig_functions")

test_that("check config", {
  config <- AnalysisConfiguration$new()
  config$fileName <- "Replicate.Name"
  config$hierarchy[["protein_Id"]] <- "Protein.Name"
  config$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  config$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
  config$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence", "Precursor.Charge", "Fragment.Ion", "Product.Charge")
  config$ident_qValue <- "annotation_QValue"
  config$set_response("Area")
  config$isotopeLabel <- "Isotope.Label"
  config$factors[["Time"]] = "Sampling.Time.Point"
  expect_equal(config$factor_keys(), "Time")
  expect_equal(config$hierarchy_keys(), c("protein_Id", "peptide_Id", "precursor_Id", "fragment_Id"))
})
