# setup raw (pre-setup_analysis) peptide data + matching config
.dup_setup_fixture <- function() {
  set.seed(1234)
  data <- prolfqua::sim_lfq_data(Nprot = 5, N = 3, PEPTIDE = TRUE)
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
  expect_error(prolfqua::setup_analysis(dup, fx$config), "more than ONE")
})

test_that("setup_analysis(debug = TRUE) returns the count table for inspection", {
  fx <- .dup_setup_fixture()
  dup <- dplyr::bind_rows(fx$data, fx$data[1, , drop = FALSE])

  expect_warning(res <- prolfqua::setup_analysis(dup, fx$config, debug = TRUE), "more than ONE")
  expect_true("n" %in% colnames(res))
  expect_gt(max(res$n), 1)
})
