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
