# setup raw (pre-setup_analysis) peptide data + matching config
.dup_setup_fixture <- function() {
  withr::local_seed(1234)
  data <- prolfqua::sim_lfq_data(n_proteins = 5, n_per_group = 3, peptide = TRUE)
  data$nr_children <- 1
  data$isotopeLabel <- "light"
  data$qValue <- 0

  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "sample"
  config$factors[["group_"]] <- "group"
  config$hierarchy[["protein_Id"]] <- c("proteinID", "idtype2")
  config$hierarchy[["peptide_Id"]] <- "peptideID"
  config$set_response("abundance")
  list(data = data, config = config)
}

test_that("setup_analysis stops on duplicate observations per key", {
  fx <- .dup_setup_fixture()

  # clean input sets up without error and has no count column
  clean <- prolfqua::setup_analysis(fx$data, fx$config)
  expect_false("n" %in% colnames(clean))

  # one duplicated observation for a (sample, protein, peptide) key
  dup <- dplyr::bind_rows(fx$data, fx$data[1, , drop = FALSE])

  # before the fix this warned and returned the count table instead of erroring
  expect_error(prolfqua::setup_analysis(dup, fx$config), "more than one observation")
})

test_that("setup_analysis(debug = TRUE) returns the count table for inspection", {
  fx <- .dup_setup_fixture()
  dup <- dplyr::bind_rows(fx$data, fx$data[1, , drop = FALSE])

  expect_warning(res <- prolfqua::setup_analysis(dup, fx$config, debug = TRUE), "more than one observation")
  expect_true("n" %in% colnames(res))
  expect_gt(max(res$n), 1)
})

test_that("setup_analysis derives sample names and default metadata", {
  fx <- .dup_setup_fixture()
  fx$data$isotopeLabel <- NULL
  fx$data$qValue <- NULL
  fx$data$nr_children <- NULL

  # isotopeLabel and qValue are optional and only reported; missing peptide
  # counts change model weights, so nr_children still warns
  expect_warning(
    messages <- capture_messages(result <- prolfqua::setup_analysis(fx$data, fx$config)),
    "no nr_children"
  )
  expect_true(any(grepl("no isotopeLabel", messages)))
  expect_true(any(grepl("no qValue", messages)))

  expect_equal(
    dplyr::n_distinct(result[[fx$config$sample_name]]),
    dplyr::n_distinct(result[[fx$config$file_name]])
  )
  expect_setequal(result[[fx$config$isotope_label]], "light")
  expect_setequal(result[[fx$config$ident_q_value]], 0)
  expect_setequal(result[[fx$config$nr_children]], 1)
})

test_that("setup_analysis works when sample_name and file_name are the same column", {
  fx <- .dup_setup_fixture()
  fx$config$sample_name <- fx$config$file_name

  result <- suppressMessages(prolfqua::setup_analysis(fx$data, fx$config))
  expect_equal(sum(colnames(result) == fx$config$file_name), 1)
  expect_equal(nrow(result), nrow(fx$data))
})

test_that("sim_lfq_data respects the mean_prot argument", {
  withr::local_seed(42)
  low <- prolfqua::sim_lfq_data(n_proteins = 50, n_per_group = 3, mean_prot = 20)
  withr::local_seed(42)
  high <- prolfqua::sim_lfq_data(n_proteins = 50, n_per_group = 3, mean_prot = 200)
  # higher mean_prot must shift abundances up; previously log(20) was hardcoded
  expect_gt(mean(high$abundance, na.rm = TRUE), mean(low$abundance, na.rm = TRUE))
})

test_that("separate_hierarchy splits hierarchy keys at the separator into their configured columns", {
  bb <- sim_lfq_data_peptide_config(n_proteins = 3)
  lfq <- LFQData$new(bb$data, bb$config)
  config <- lfq$get_config()$clone(deep = TRUE)
  data <- lfq$data_long()
  data$protein_Id[1] <- "noseparator"

  res <- separate_hierarchy(data, config)
  expect_identical(colnames(res)[5:7], c("protein_Id", "proteinID", "idtype2"))
  parts <- strsplit(data$protein_Id[-1], "~", fixed = TRUE)
  expect_identical(res$proteinID[-1], vapply(parts, `[`, character(1), 1))
  expect_identical(res$idtype2[-1], vapply(parts, `[`, character(1), 2))
  # too few parts fill from the start
  expect_identical(res$proteinID[1], "noseparator")
  expect_true(is.na(res$idtype2[1]))
  # keys below the hierarchy depth are left alone
  expect_false("peptideID" %in% colnames(res))

  config$hierarchy_depth <- 2
  deeper <- separate_hierarchy(data, config)
  expect_identical(deeper$peptideID, data$peptide_Id)
})
