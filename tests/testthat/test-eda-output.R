# Default print and plot methods, numeric precision of printed output (srr EDA standards).

eda_fixture <- function() {
  sim <- sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
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

test_that("default plots use the prolfqua theme, with gridlines only on fitted-value plots", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- eda_fixture()
  protein <- suppressWarnings(lfq$get_Aggregator()$aggregate())
  fit <- build_contrast_analysis(protein, "~ group_", c(AvsC = "group_A - group_Ctrl"), "lm")
  themed <- function(p) {
    th <- p$theme
    expect_identical(ggplot2::calc_element("panel.background", th)$fill, "white")
    expect_identical(ggplot2::calc_element("axis.line.x", th)$colour, "black")
    expect_identical(ggplot2::calc_element("axis.title.x", th)$face, "bold")
    inherits(ggplot2::calc_element("panel.grid.major.x", th), "element_line")
  }
  expect_false(themed(plot(lfq)))
  expect_true(themed(plot(fit)))
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

test_that("duplicate sample and hierarchy keys are rejected at data, model and aggregation entry", {
  lfq <- eda_fixture()
  data <- lfq$data_long()
  expect_error(
    LFQData$new(rbind(data, data[1:3, ]), lfq$get_config()),
    "More than one row for 3",
    class = "prolfqua_error_duplicate_keys"
  )
  expect_error(LFQData$new(data[, -1], lfq$get_config()), class = "prolfqua_error_missing_column")
  doubled <- lfq$get_copy()$set_data(rbind(data, data[1:3, ]))
  expect_error(
    build_contrast_analysis(doubled, "~ group_ + (1 | peptide_Id)", c(AvsC = "group_A - group_Ctrl"), "lmer_nested"),
    class = "prolfqua_error_duplicate_keys"
  )
  expect_error(doubled$get_Aggregator()$aggregate(), "More than one row for 3", class = "prolfqua_error_duplicate_keys")
})

test_that("a response with extra attributes is used as plain numeric values", {
  withr::local_options(prolfqua.verbose = TRUE)
  sim <- sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
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
  sim <- sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
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

eda_tables <- function() {
  withr::local_options(prolfqua.verbose = FALSE)
  peptide <- eda_fixture()
  protein <- suppressWarnings(peptide$get_Aggregator()$aggregate())
  fit <- build_contrast_analysis(
    protein,
    "~ group_",
    c(AvsC = "group_A - group_Ctrl", BvsC = "group_B - group_Ctrl"),
    "lm"
  )
  list(peptide = peptide, protein = protein, fit = fit)
}

column_classes <- function(table) vapply(table, function(column) class(column)[1], character(1))

#' @srrstats {EA6.0, EA6.0b, EA6.0d} Main returned tables (contrasts, coefficients, per-group statistics, variance
#'   explained, hierarchy counts, provenance, sample annotation) are tested for dimensions and the class of every
#'   column; with the class, column-name and value tests (EA6.0a, EA6.0c, EA6.0e) this covers EA6.0.
#' @noRd
NULL

test_that("main returned tables have the documented dimensions and column classes", {
  fixture <- eda_tables()
  n_protein <- nrow(fixture$protein$hierarchy())
  n_peptide <- nrow(unique(fixture$peptide$data_long()[, fixture$peptide$hierarchy_keys()]))

  contrasts <- fixture$fit$get_contrasts()
  expect_identical(dim(contrasts), c(2L * n_protein, 16L))
  expect_identical(
    unname(column_classes(contrasts)),
    c(rep("character", 4), "numeric", "numeric", "numeric", "integer", rep("numeric", 8))
  )

  coefs <- coef(fixture$fit)
  expect_identical(dim(coefs), c(3L * n_protein, 5L))
  expect_identical(unname(column_classes(coefs)), c("character", "character", "numeric", "numeric", "numeric"))

  stats <- fixture$peptide$get_Stats()$stats()
  expect_identical(dim(stats), c(4L * n_peptide, 11L)) # three groups and "All"
  expect_identical(
    unname(column_classes(stats)),
    c(rep("character", 3), rep("integer", 3), rep("numeric", 4), "character")
  )

  variance <- summarize_variance_explained(fixture$protein)
  expect_identical(dim(variance), c(2L * n_protein, 7L))
  expect_identical(unname(column_classes(variance)), c(rep("character", 2), "integer", rep("numeric", 4)))

  counts <- fixture$peptide$hierarchy_counts()
  expect_identical(dim(counts), c(1L, 2L))
  expect_identical(unname(column_classes(counts)), c("integer", "integer"))
  expect_identical(counts$protein_Id, n_protein)

  provenance <- fixture$protein$provenance()
  expect_identical(ncol(provenance), 7L)
  expect_identical(
    unname(column_classes(provenance)),
    c("character", "character", "integer", "list", "character", "logical", "character")
  )

  factors <- fixture$peptide$factors()
  expect_identical(dim(factors), c(12L, 3L))
  expect_true(all(column_classes(factors) == "character"))
})

#' @srrstats {EA6.1} Default plots are tested for their layers, aesthetic mappings, axis titles and colour scales by
#'   building them with ggplot2 (no image snapshots: vdiffr is not used, so tests stay stable across ggplot2
#'   versions).
#' @noRd
NULL

test_that("default plots have the expected layers, mappings, labels and colour scales", {
  fixture <- eda_tables()
  layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
  # missing intensities are dropped from the distribution plots with ggplot2's "Removed ... rows" warning
  build <- function(p) {
    withCallingHandlers(
      ggplot2::ggplot_build(p),
      warning = function(w) if (grepl("^Removed [0-9]+ rows", conditionMessage(w))) invokeRestart("muffleWarning")
    )
  }

  violin <- plot(fixture$peptide)
  build(violin)
  expect_true("GeomViolin" %in% layer_geoms(violin))
  expect_identical(rlang::as_label(violin$mapping$y), fixture$peptide$response())

  plotter <- fixture$peptide$get_Plotter()
  density <- plotter$intensity_distribution_density()
  build(density)
  expect_identical(rlang::as_label(density$mapping$x), fixture$peptide$response())

  pca <- fixture$protein$get_Plotter()$pca()
  build(pca)
  expect_match(pca$labels$x %||% ggplot2::get_labs(pca)$x, "^PC1 \\([0-9]+% variance\\)$")
  expect_true("GeomPoint" %in% layer_geoms(pca))

  volcano <- plot(fixture$fit, type = "volcano")
  built <- build(volcano)
  expect_true("GeomPoint" %in% layer_geoms(volcano))
  expect_identical(
    built$plot$scales$get_scales("colour")$map("observed"),
    "#000000"
  )
  expect_true(all(c("GeomHline", "GeomVline") %in% layer_geoms(volcano)))

  fitted <- plot(fixture$fit)
  build(fitted)
  expect_s3_class(fitted, "ggplot")
  expect_gt(length(fitted$layers), 1L)
})
