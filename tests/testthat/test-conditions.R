# Typed error conditions raised at prolfqua's public boundaries.

test_that("setup_analysis aborts (typed) when the file_name column is missing", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 5, seed = 42)
  cfg <- istar$config$clone(deep = TRUE)
  cfg$file_name <- "column_that_does_not_exist"

  expect_error(
    prolfqua::setup_analysis(istar$data, cfg),
    class = "prolfqua_error_missing_column"
  )
})

test_that("setup_analysis aborts (typed) when no factors are configured", {
  raw <- data.frame(
    sample = c("s1", "s2"),
    proteinID = c("P1", "P1"),
    idtype2 = c("x", "x"),
    peptideID = c("pep1", "pep1"),
    abundance = c(10, 12),
    stringsAsFactors = FALSE
  )
  cfg <- prolfqua::AnalysisConfiguration$new()
  cfg$file_name <- "sample"
  cfg$hierarchy[["protein_Id"]] <- c("proteinID", "idtype2")
  cfg$hierarchy[["peptide_Id"]] <- "peptideID"
  cfg$set_response("abundance")

  expect_error(
    prolfqua::setup_analysis(raw, cfg),
    class = "prolfqua_error_invalid_configuration"
  )
})

test_that("setup_analysis aborts (typed) on duplicate hierarchy-key / sample combinations", {
  raw <- data.frame(
    sample = c("s1", "s1"),
    group = c("A", "A"),
    proteinID = c("P1", "P1"),
    idtype2 = c("x", "x"),
    peptideID = c("pep1", "pep1"),
    abundance = c(10, 12),
    isotopeLabel = "light",
    qValue = 0,
    stringsAsFactors = FALSE
  )
  cfg <- prolfqua::AnalysisConfiguration$new()
  cfg$file_name <- "sample"
  cfg$factors["group_"] <- "group"
  cfg$hierarchy[["protein_Id"]] <- c("proteinID", "idtype2")
  cfg$hierarchy[["peptide_Id"]] <- "peptideID"
  cfg$set_response("abundance")

  # raw has no nr_children column, which setup_analysis() warns about before it checks the keys
  expect_error(
    muffle_warnings_matching(prolfqua::setup_analysis(raw, cfg), "no nr_children column"),
    class = "prolfqua_error_duplicate_keys"
  )
})

test_that("LFQData$set_data rejects a non-data-frame (typed)", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_error(
    lfq$set_data(list(a = 1)),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("LFQData$set_data enforces the configured schema (typed)", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_error(
    lfq$set_data(data.frame(unrelated = 1)),
    class = "prolfqua_error_missing_column"
  )
})

test_that("build_contrast_analysis aborts (typed) on an unknown method", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  lfq$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  expect_error(
    prolfqua::build_contrast_analysis(lfq, "~ group_", contrasts, method = "definitely_not_real"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("LFQData$get_Aggregator aborts (typed) on an unknown method", {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_error(
    lfq$get_Aggregator("not_a_method"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("aggregating single-level hierarchy aborts (typed)", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_error(
    lfq$get_Aggregator("medpolish"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("center_to_reference aborts (typed) on mismatched transformation state", {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  ref <- lfq$get_copy()
  ref <- ref$get_Transformer()$log2()$lfq

  expect_error(
    lfq$get_Transformer()$center_to_reference(ref),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("informative messages carry class prolfqua_message", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_message(
    lfq$get_Transformer()$log2(),
    "Column added",
    class = "prolfqua_message"
  )
})

test_that("options(prolfqua.verbose = FALSE) silences informative messages", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  withr::local_options(prolfqua.verbose = FALSE)

  expect_no_message(lfq$get_Transformer()$log2())
})

test_that("errors carry class prolfqua_error", {
  expect_error(
    prolfqua::register_facade(name = 1, class = "X", needs = "aggregated"),
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    prolfqua::contrasts_linfct(data.frame(protein_Id = c("p1", "p2")), linfct = list(1)),
    class = "prolfqua_error"
  )
})

test_that("plot_pca rejects a PC argument of the wrong length", {
  x <- matrix(rnorm(20), nrow = 5, dimnames = list(NULL, paste0("s", 1:4)))
  annot <- data.frame(sample = paste0("s", 1:4), group = c("a", "a", "b", "b"))
  expect_error(
    prolfqua::plot_pca(x, annot, sample_name = "sample", factor_keys = "group", PC = 1),
    class = "prolfqua_error_bad_argument"
  )
})
