# Test ContrastsPlotter R6 class

test_that("ContrastsPlotter produces correct plot types", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20)
  config <- istar$config

  modelFunction <- strategy_lmer(
    "abundance ~ group_ + (1 | peptide_Id) + (1 | sample)"
  )
  mod <- build_model(istar$data, modelFunction, subject_id = config$hierarchy_keys_depth())

  Contr <- c(
    "groupA_vs_Ctrl" = "group_A - group_Ctrl",
    "groupB_vs_Ctrl" = "group_B - group_Ctrl"
  )
  contrast <- Contrasts$new(mod, Contr)
  contr <- contrast$get_contrasts()

  cp <- ContrastsPlotter$new(
    contr,
    contrast$subject_id,
    volcano = list(list(score = "FDR", thresh = 0.1)),
    histogram = list(
      list(score = "p.value", xlim = c(0, 1, 0.05)),
      list(score = "FDR", xlim = c(0, 1, 0.05))
    ),
    score = list(list(score = "statistic", thresh = 5))
  )
  expect_setequal(intersect("contrast_df", names(cp)), character())

  # volcano returns list of ggplots
  res <- cp$volcano()
  expect_type(res, "list")
  expect_s3_class(res$FDR, "ggplot")

  # volcano_plotly returns list of plotly
  res_plotly <- cp$volcano_plotly()
  expect_type(res_plotly, "list")
  expect_s3_class(res_plotly$FDR, "plotly")

  # ma_plot returns ggplot
  ma <- cp$ma_plot()
  expect_s3_class(ma, "ggplot")

  # ma_plotly returns plotly
  ma_pl <- cp$ma_plotly(rank = TRUE)
  expect_s3_class(ma_pl, "plotly")

  # histogram returns list of ggplots
  h <- cp$histogram()
  expect_type(h, "list")
  expect_s3_class(h$FDR, "ggplot")
  expect_s3_class(h$p.value, "ggplot")

  # histogram_estimate returns ggplot
  he <- cp$histogram_estimate()
  expect_s3_class(he, "ggplot")

  # histogram_diff returns ggplot
  hd <- cp$histogram_diff()
  expect_s3_class(hd, "ggplot")

  # score_plot returns list of ggplots
  sp <- cp$score_plot(legend = FALSE)
  expect_type(sp, "list")
  expect_s3_class(sp$statistic, "ggplot")

  # barplot_threshold returns list with plot element
  bt <- cp$barplot_threshold()
  expect_type(bt, "list")
  expect_s3_class(bt$FDR$plot, "ggplot")
})

test_that("volcano score names can differ from score columns", {
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
