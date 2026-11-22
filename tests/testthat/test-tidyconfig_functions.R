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
  config$factors[["Time"]] <- "Sampling.Time.Point"
  expect_equal(config$factor_keys(), "Time")
  expect_equal(config$hierarchy_keys(), c("protein_Id", "peptide_Id", "precursor_Id", "fragment_Id"))
})

test_that("pattern_decoys / pattern_contaminants default to NULL and round-trip", {
  config <- AnalysisConfiguration$new()
  expect_null(config$pattern_decoys)
  expect_null(config$pattern_contaminants)

  # NULL (unset) round-trip: absent from extracted list, restored as default NULL
  vals <- r6_extract_values(config)
  expect_false("pattern_decoys" %in% names(vals))
  rebuilt <- list_to_analysis_configuration(vals)
  expect_null(rebuilt$pattern_decoys)
  expect_null(rebuilt$pattern_contaminants)

  # configured patterns survive an extract -> reconstruct round-trip
  config$pattern_decoys <- "^shuffled_"
  config$pattern_contaminants <- "^KERATIN_"
  vals2 <- r6_extract_values(config)
  expect_identical(vals2$pattern_decoys, "^shuffled_")
  expect_identical(vals2$pattern_contaminants, "^KERATIN_")
  rebuilt2 <- list_to_analysis_configuration(vals2)
  expect_identical(rebuilt2$pattern_decoys, "^shuffled_")
  expect_identical(rebuilt2$pattern_contaminants, "^KERATIN_")
})

test_that("list_to_contrast_configuration restores roles and methods", {
  cfg <- ContrastConfiguration$new(
    subject_id = "protein_Id",
    contrast_col = "Bait",
    effect_col = "log2_EFCs",
    score_col = "SaintScore",
    pvalue_col = NA_character_,
    fdr_col = "BFDR",
    significance_directional = TRUE
  )

  vals <- r6_extract_values(cfg)
  rebuilt <- list_to_contrast_configuration(vals)

  expect_identical(r6_extract_values(rebuilt), vals)
  # r6_extract_values drops the methods; the reconstruction brings them back.
  expect_false(rebuilt$has_pvalue())
  expect_true(rebuilt$significance_directional)
  expect_true(ContrastConfiguration$new()$has_pvalue())
  # An unset field keeps the schema default rather than becoming NULL.
  expect_identical(
    list_to_contrast_configuration(list())$fdr_col,
    ContrastConfiguration$new()$fdr_col
  )
  expect_error(
    list_to_contrast_configuration(list(no_such_field = "x")),
    "Unknown ContrastConfiguration fields"
  )
})

test_that("make_interaction_column keeps its levels and values (review B5)", {
  xx <- data.frame(A = c("a", "b", "a"), B = c("d", "d", "e"))
  expect_identical(
    make_interaction_column(xx, c("B", "A"))$interaction,
    factor(c("Bd.Aa", "Bd.Ab", "Be.Aa"), levels = c("Bd.Aa", "Be.Aa", "Bd.Ab", "Be.Ab"))
  )
  expect_identical(make_interaction_column(xx, "A")$interaction, factor(c("Aa", "Ab", "Aa")))
})

test_that("table_facade returns a kable table, or the data frame with kable = FALSE", {
  df <- head(iris, 3)
  expect_identical(table_facade(df, caption = "x", kable = FALSE), df)
  skip_if_not_installed("knitr")
  expect_s3_class(table_facade(df, caption = "x"), "knitr_kable")
})

test_that("list_to_analysis_configuration names unknown keys and method names in a classed error", {
  vals <- r6_extract_values(AnalysisConfiguration$new())
  expect_error(
    list_to_analysis_configuration(c(vals, list(no_such_field = 1))),
    "no_such_field",
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    list_to_analysis_configuration(list(hierarchy_keys = "x")),
    "hierarchy_keys",
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    list_to_analysis_configuration(list(table = list(not_a_field = 1))),
    "not_a_field",
    class = "prolfqua_error_bad_argument"
  )
  # the legacy parameter block keeps only known fields
  legacy <- list_to_analysis_configuration(list(table = list(sep = "~"), parameter = list(old_parameter = 1)))
  expect_identical(legacy$sep, "~")
})

test_that("list_to_analysis_configuration loads configurations stored with the removed isotope_label", {
  vals <- r6_extract_values(AnalysisConfiguration$new())
  expect_false("isotope_label" %in% names(vals))
  flat <- list_to_analysis_configuration(c(vals, list(isotope_label = "isotopeLabel")))
  expect_identical(r6_extract_values(flat), vals)
  nested <- list_to_analysis_configuration(list(table = list(sep = "~", isotope_label = "isotopeLabel")))
  expect_identical(nested$sep, "~")
})
