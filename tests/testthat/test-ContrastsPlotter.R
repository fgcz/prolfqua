# Test ContrastsPlotter R6 class

# lm contrasts of 20 simulated proteins for two contrasts
plotter_contrasts <- function() {
  lfqdata <- protein_lfqdata_fixture(n_proteins = 20)
  contrasts <- c(
    "groupA_vs_Ctrl" = "group_A - group_Ctrl",
    "groupB_vs_Ctrl" = "group_B - group_Ctrl"
  )
  build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")$get_contrasts()
}

new_contrasts_plotter <- function(contr = plotter_contrasts()) {
  ContrastsPlotter$new(
    contr,
    "protein_Id",
    volcano = list(list(score = "FDR", thresh = 0.1)),
    histogram = list(
      list(score = "p.value", breaks = seq(0, 1, by = 0.05)),
      list(score = "FDR", breaks = seq(0, 1, by = 0.1))
    ),
    score = list(list(score = "statistic", thresh = 5))
  )
}

test_that("ContrastsPlotter does not expose the contrast table as a field", {
  cp <- new_contrasts_plotter()
  expect_setequal(intersect("contrast_df", names(cp)), character())
})

test_that("ContrastsPlotter volcano returns a list of ggplots", {
  res <- new_contrasts_plotter()$volcano()
  expect_type(res, "list")
  expect_s3_class(res$FDR, "ggplot")
})

test_that("ContrastsPlotter volcano_plotly returns a list of plotly objects", {
  skip_if_not_installed("plotly")
  res_plotly <- new_contrasts_plotter()$volcano_plotly()
  expect_type(res_plotly, "list")
  expect_s3_class(res_plotly$FDR, "plotly")
})

test_that("ContrastsPlotter ma_plot returns a ggplot", {
  expect_s3_class(new_contrasts_plotter()$ma_plot(), "ggplot")
})

test_that("ContrastsPlotter ma_plotly returns a plotly object", {
  skip_if_not_installed("plotly")
  expect_s3_class(new_contrasts_plotter()$ma_plotly(rank = TRUE), "plotly")
})

test_that("ContrastsPlotter histogram returns a list of ggplots with the configured breaks", {
  contr <- plotter_contrasts()
  h <- new_contrasts_plotter(contr)$histogram()
  expect_type(h, "list")
  expect_s3_class(h$FDR, "ggplot")
  expect_s3_class(h$p.value, "ggplot")
  expect_equal(h$FDR$layers[[1]]$stat_params$breaks, seq(0, 1, by = 0.1))
  # the old (from, to, by) triple under xlim is no longer read
  old_spec <- ContrastsPlotter$new(
    contr,
    "protein_Id",
    histogram = list(list(score = "FDR", xlim = c(0, 1, 0.05)))
  )
  expect_error(old_spec$histogram(), class = "prolfqua_error_bad_argument")
})

test_that("ContrastsPlotter histogram_diff returns a ggplot", {
  expect_s3_class(new_contrasts_plotter()$histogram_diff(), "ggplot")
})

test_that("ContrastsPlotter score_plot returns a list of ggplots", {
  sp <- new_contrasts_plotter()$score_plot(legend = FALSE)
  expect_type(sp, "list")
  expect_s3_class(sp$statistic, "ggplot")
})

test_that("ContrastsPlotter barplot_threshold returns a list with a plot element", {
  bt <- new_contrasts_plotter()$barplot_threshold()
  expect_type(bt, "list")
  expect_s3_class(bt$FDR$plot, "ggplot")
})

test_that("volcano score names can differ from score columns", {
  skip_if_not_installed("plotly")
  contrast_df <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    contrast = "A_vs_B",
    modelName = "model",
    diff = c(-2, 0, 2),
    BFDR = c(0.01, 0.1, 0.2)
  )
  cp <- ContrastsPlotter$new(
    contrast_df,
    subject_id = "protein_Id",
    volcano = list(list(score = "BFDR", name = "FDR", thresh = 0.1))
  )

  volcano <- cp$volcano()
  expect_named(volcano, "FDR")
  expect_s3_class(volcano$FDR, "ggplot")
  expect_equal(volcano$FDR$labels$y, "-log10(BFDR)")

  volcano_plotly <- cp$volcano_plotly()
  expect_named(volcano_plotly, "FDR")
  expect_s3_class(volcano_plotly$FDR, "plotly")
})

test_that("estimate types use stable colours across contrast plots", {
  contrast_df <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    contrast = "A_vs_B",
    modelName = "lm_impute",
    estimate_type = c("observed", "lod_imputed", "missing_fallback"),
    diff = c(-2, 0, 2),
    FDR = c(0.01, 0.1, 0.2),
    statistic = c(-3, 0, 3),
    avgAbd = c(10, 11, 12)
  )
  cp <- ContrastsPlotter$new(
    contrast_df,
    subject_id = "protein_Id",
    volcano = list(list(score = "FDR", thresh = 0.1)),
    score = list(list(score = "statistic", thresh = 2))
  )
  plots <- list(
    cp$volcano()$FDR,
    cp$ma_plot(rank = FALSE),
    cp$score_plot()$statistic
  )

  expected <- c(observed = "#000000", lod_imputed = "#E69F00", missing_fallback = "#0072B2")
  for (plot in plots) {
    built <- ggplot2::ggplot_build(plot)
    colour_scale <- built$plot$scales$get_scales("colour")
    expect_equal(colour_scale$map(names(expected)), unname(expected))
  }
})

test_that("volcano_plotly shares the crosstalk group with linked tables", {
  skip_if_not_installed("plotly")
  skip_if_not_installed("crosstalk")
  # Regression: the DT results table uses crosstalk::SharedData group "BB" so
  # brushing the volcano highlights table rows (and vice versa). The plotly
  # volcano must publish into the same group, keyed on the united subject_id.
  contrast_df <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    contrast = "A_vs_B",
    modelName = "model",
    diff = c(-2, 0, 2),
    FDR = c(0.01, 0.1, 0.2)
  )
  cp <- ContrastsPlotter$new(
    contrast_df,
    subject_id = "protein_Id",
    volcano = list(list(score = "FDR", thresh = 0.1))
  )
  expect_equal(cp$group, "BB")

  built <- plotly::plotly_build(cp$volcano_plotly()$FDR)
  expect_equal(built$x$data[[1]]$set, "BB")
  expect_setequal(unlist(built$x$data[[1]]$key), c("P1", "P2", "P3"))

  # group is configurable for callers that link against a different set
  cp2 <- ContrastsPlotter$new(
    contrast_df,
    subject_id = "protein_Id",
    volcano = list(list(score = "FDR", thresh = 0.1)),
    group = "OTHER"
  )
  built2 <- plotly::plotly_build(cp2$volcano_plotly()$FDR)
  expect_equal(built2$x$data[[1]]$set, "OTHER")
})

test_that("volcano plots do not cap small non-zero FDR values by default", {
  contrast_df <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    contrast = "A_vs_B",
    modelName = "model",
    diff = c(-2, 0, 2),
    FDR = c(1e-8, 1e-3, 0)
  )
  cp <- ContrastsPlotter$new(
    contrast_df,
    subject_id = "protein_Id",
    volcano = list(list(score = "FDR", thresh = 0.1))
  )

  default_plot <- cp$volcano()$FDR
  default_y <- ggplot2::ggplot_build(default_plot)$data[[1]]$y

  expect_gt(max(default_y[is.finite(default_y)]), 7.9)
  expect_true(all(is.finite(default_y)))

  capped_plot <- cp$volcano(min_score = 1e-4)$FDR
  capped_y <- ggplot2::ggplot_build(capped_plot)$data[[1]]$y

  expect_equal(max(capped_y), 4)
})
