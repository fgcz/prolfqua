context("test-tidyconfig_functions")

test_that("check config", {
  config <- AnalysisConfiguration$new()
  config$file_name <- "Replicate.Name"
  config$hierarchy[["protein_Id"]] <- "Protein.Name"
  config$hierarchy[["peptide_Id"]] <- "Peptide.Sequence"
  config$hierarchy[["precursor_Id"]] <- c("Peptide.Sequence", "Precursor.Charge")
  config$hierarchy[["fragment_Id"]] <- c("Peptide.Sequence", "Precursor.Charge", "Fragment.Ion", "Product.Charge")
  config$ident_q_value <- "annotation_QValue"
  config$set_response("Area")
  config$isotope_label <- "Isotope.Label"
  config$factors[["Time"]] <- "Sampling.Time.Point"
  expect_equal(config$factor_keys(), "Time")
  expect_equal(config$hierarchy_keys(), c("protein_Id", "peptide_Id", "precursor_Id", "fragment_Id"))
})

test_that("pattern_decoys / pattern_contaminants default to NULL and round-trip", {
  config <- AnalysisConfiguration$new()
  expect_null(config$pattern_decoys)
  expect_null(config$pattern_contaminants)

  # NULL (unset) round-trip: absent from extracted list, restored as default NULL
  vals <- R6_extract_values(config)
  expect_false("pattern_decoys" %in% names(vals))
  rebuilt <- list_to_AnalysisConfiguration(vals)
  expect_null(rebuilt$pattern_decoys)
  expect_null(rebuilt$pattern_contaminants)

  # configured patterns survive an extract -> reconstruct round-trip
  config$pattern_decoys <- "^shuffled_"
  config$pattern_contaminants <- "^KERATIN_"
  vals2 <- R6_extract_values(config)
  expect_identical(vals2$pattern_decoys, "^shuffled_")
  expect_identical(vals2$pattern_contaminants, "^KERATIN_")
  rebuilt2 <- list_to_AnalysisConfiguration(vals2)
  expect_identical(rebuilt2$pattern_decoys, "^shuffled_")
  expect_identical(rebuilt2$pattern_contaminants, "^KERATIN_")
})

test_that("list_to_ContrastConfiguration restores roles and methods", {
  cfg <- ContrastConfiguration$new(
    subject_id = "protein_Id",
    contrast_col = "Bait",
    effect_col = "log2_EFCs",
    score_col = "SaintScore",
    pvalue_col = NA_character_,
    fdr_col = "BFDR",
    significance_directional = TRUE
  )

  vals <- R6_extract_values(cfg)
  rebuilt <- list_to_ContrastConfiguration(vals)

  expect_identical(R6_extract_values(rebuilt), vals)
  # R6_extract_values drops the methods; the reconstruction brings them back.
  expect_false(rebuilt$has_pvalue())
  expect_true(rebuilt$significance_directional)
  expect_true(ContrastConfiguration$new()$has_pvalue())
  # An unset field keeps the schema default rather than becoming NULL.
  expect_identical(
    list_to_ContrastConfiguration(list())$fdr_col,
    ContrastConfiguration$new()$fdr_col
  )
  expect_error(
    list_to_ContrastConfiguration(list(no_such_field = "x")),
    "Not ContrastConfiguration fields"
  )
})
