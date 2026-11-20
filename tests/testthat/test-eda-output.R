# Default print and plot methods, numeric precision of printed output (srr EDA standards).

eda_fixture <- function() {
  sim <- sim_lfq_data_peptide_config(Nprot = 5, seed = 42)
  suppressMessages(LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq)
}

test_that("LFQData prints samples, response, hierarchy, factor levels, missing values and provenance", {
  lfq <- eda_fixture()
  out <- capture.output(print(lfq))
  expect_match(out[1], "^LFQData: 12 samples$")
  expect_match(out[2], "log2_abundance \\(transformed\\)")
  expect_match(out[3], "protein_Id 5, peptide_Id")
  expect_match(out[4], "^Samples per group_: A 4, B 4, Ctrl 4$")
  expect_match(out[5], "^Response range: .*; missing [0-9]+ of [0-9]+ \\(.*%\\)$")
  expect_match(out[6], "^Provenance:")
  expect_invisible(print(lfq))
  range_text <- function(digits) {
    sub("^Response range: ([^;]+);.*$", "\\1", capture.output(print(lfq, digits = digits))[5])
  }
  expect_match(range_text(2), "^[0-9.]{1,3} to [0-9.]{1,3}$")
  expect_false(identical(range_text(2), range_text(6)))
})

test_that("plot() of LFQData shows the response distribution per sample", {
  lfq <- eda_fixture()
  p <- plot(lfq)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$layers[[1]]$geom, "GeomViolin"))
})

#' @srrstats {EA5.3} Summary tables are tibbles, whose printed header states the class of each column (`<int>`,
#'   `<dbl>`, `<chr>`); tested for the analysis summary.
#' @noRd
NULL

test_that("analysis print and summary round to digits and show column classes", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- suppressWarnings(eda_fixture()$get_Aggregator()$aggregate())
  fit <- build_contrast_analysis(protein, "~ group_", c(AvsC = "group_A - group_Ctrl"), "lm")
  coef_lines <- function(digits) {
    out <- capture.output(print(fit, digits = digits))
    out[grep("(Intercept)", out, fixed = TRUE)[1]]
  }
  short <- regmatches(coef_lines(2), gregexpr("[0-9]+\\.[0-9]+", coef_lines(2)))[[1]]
  expect_true(all(nchar(gsub("^0\\.0*|\\.", "", short)) <= 2))
  expect_false(identical(coef_lines(2), coef_lines(6)))
  summary_out <- capture.output(print(summary(fit), digits = 2))
  expect_true(any(grepl("<int>", summary_out)) && any(grepl("<dbl>", summary_out)))
  expect_s3_class(summary(fit)$per_contrast, "tbl_df")
  table3 <- fit$contrast$contrast_summary_table()
  table6 <- fit$contrast$contrast_summary_table(digits = 6)
  expect_equal(table3$effect, signif(table6$effect, 3))
  expect_false(isTRUE(all.equal(table3$effect, table6$effect, tolerance = 0)))
})

test_that("duplicate sample and hierarchy keys are rejected at model and aggregation entry", {
  lfq <- eda_fixture()
  data <- lfq$data_long()
  doubled <- LFQData$new(rbind(data, data[1:3, ]), lfq$get_config())
  expect_error(
    build_contrast_analysis(doubled, "~ group_ + (1 | peptide_Id)", c(AvsC = "group_A - group_Ctrl"), "lmer_nested"),
    class = "prolfqua_error_duplicate_keys"
  )
  expect_error(doubled$get_Aggregator()$aggregate(), "More than one row for 3", class = "prolfqua_error_duplicate_keys")
})

test_that("a response with extra attributes is used as plain numeric values", {
  withr::local_options(prolfqua.verbose = TRUE)
  sim <- sim_lfq_data_peptide_config(Nprot = 5, seed = 42)
  plain <- LFQData$new(sim$data, sim$config)
  for (wrap in list(I, function(x) structure(x, label = "Intensity", units = "counts"))) {
    data <- sim$data
    data$abundance <- wrap(data$abundance)
    expect_message(wrapped <- LFQData$new(data, sim$config), "dropped attributes", class = "prolfqua_message")
    expect_null(attributes(wrapped$data_long()$abundance))
    expect_identical(wrapped$data_long()$abundance, plain$data_long()$abundance)
    log2_values <- suppressMessages(wrapped$get_Transformer()$log2()$lfq)$data_long()$log2_abundance
    expect_type(log2_values, "double")
    expect_null(attributes(log2_values))
  }
})

#' @srrstats {EA4.0} Tables come back as data frames (tibbles where a step computes a new table) whatever data frame
#'   class went in; numeric responses stay double through transformation, aggregation and contrasts.
#' @noRd
NULL

test_that("return types follow the input: data frames in, data frames and double responses out", {
  withr::local_options(prolfqua.verbose = FALSE)
  sim <- sim_lfq_data_peptide_config(Nprot = 5, seed = 42)
  lfq <- LFQData$new(as.data.frame(sim$data), sim$config)
  expect_s3_class(lfq$data_long(), "data.frame")
  log2 <- lfq$get_Transformer()$log2()$lfq
  protein <- suppressWarnings(log2$get_Aggregator()$aggregate())
  fit <- build_contrast_analysis(protein, "~ group_", c(AvsC = "group_A - group_Ctrl"), "lm")
  tables <- list(log2$data_long(), protein$data_long(), protein$data_wide()$data, fit$get_contrasts(), log2$factors())
  for (table in tables) {
    expect_s3_class(table, "data.frame")
  }
  expect_type(log2$data_long()[[log2$response()]], "double")
  expect_type(protein$data_long()[[protein$response()]], "double")
  expect_type(fit$get_contrasts()$diff, "double")
})
