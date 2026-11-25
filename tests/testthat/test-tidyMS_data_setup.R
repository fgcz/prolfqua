# setup raw (pre-setup_analysis) peptide data + matching config
.dup_setup_fixture <- function() {
  withr::local_seed(1234)
  data <- prolfqua::sim_lfq_data(n_proteins = 5, n_per_group = 3, peptide = TRUE)
  data$nr_children <- 1
  data$identification_confidence <- 1

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

test_that("setup_analysis rejects two isotope labels as duplicate observations", {
  fx <- .dup_setup_fixture()
  two_labels <- dplyr::bind_rows(
    dplyr::mutate(fx$data, isotopeLabel = "light"),
    dplyr::mutate(fx$data, isotopeLabel = "heavy")
  )
  expect_error(prolfqua::setup_analysis(two_labels, fx$config), "more than one observation")
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
  fx$data$identification_confidence <- NULL
  fx$data$nr_children <- NULL

  # a missing identification confidence defaults to 1 with a message; missing peptide
  # counts change model weights, so nr_children still warns
  expect_warning(
    messages <- capture_messages(result <- prolfqua::setup_analysis(fx$data, fx$config)),
    "no nr_children"
  )
  expect_true(any(grepl("no identification_confidence", messages)))
  expect_true(any(grepl("derived identification_confidence_feature, nr_children_feature", messages)))

  expect_equal(
    dplyr::n_distinct(result[[fx$config$sample_name]]),
    dplyr::n_distinct(result[[fx$config$file_name]])
  )
  expect_setequal(result[[fx$config$identification_confidence]], 1)
  expect_setequal(result[[fx$config$identification_confidence_feature]], 1)
  expect_setequal(result[[fx$config$nr_children]], 1)
  expect_setequal(result[[fx$config$nr_children_feature]], 1)
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
