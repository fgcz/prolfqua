# Plot correctness: per-panel summaries, colour scales, non-syntactic column names, point sizes.

.plot_fix_contrasts <- function() {
  istar <- sim_lfq_data_peptide_config(n_proteins = 20, with_missing = FALSE, seed = 1)
  lfq <- LFQData$new(istar$data, istar$config)$get_Transformer()$log2()$lfq
  agg <- lfq$get_Aggregator()
  agg$aggregate()
  facade <- build_contrast_analysis(
    agg$lfq_agg,
    "~ group_",
    c(A_vs_Ctrl = "group_A - group_Ctrl", B_vs_Ctrl = "group_B - group_Ctrl"),
    method = "lm"
  )
  facade$get_contrasts()
}

.volcano_input <- function(n_colours) {
  withr::with_seed(
    2,
    data.frame(
      logFC = stats::rnorm(60),
      `adj p` = stats::runif(60),
      Condition = rep(c("a", "b"), 30),
      colour = rep(letters[seq_len(n_colours)], length.out = 60),
      Name = paste0("n", 1:60),
      check.names = FALSE
    )
  )
}

test_that("histogram_diff() draws each contrast's own median once per panel", {
  contrasts <- .plot_fix_contrasts()
  cp <- ContrastsPlotter$new(contrasts, "protein_Id")
  vline <- ggplot2::layer_data(cp$histogram_diff(), 2)
  medians <- tapply(contrasts$diff, contrasts$contrast, stats::median, na.rm = TRUE)
  expect_equal(nrow(vline), length(medians))
  panels <- as.character(ggplot2::ggplot_build(cp$histogram_diff())$layout$layout$contrast)
  expect_equal(vline$xintercept[order(vline$PANEL)], as.vector(medians[panels]))
})

test_that("volcano colours: up to 4 levels keep the palette, more levels still build", {
  expect_equal(.contrast_colour_values("colour", letters[1:3]), c("#000000", "#E69F00", "#0072B2", "#D55E00"))
  expect_length(.contrast_colour_values("colour", letters[1:6]), 6)
  p <- .multigroup_volcano_labelled(.volcano_input(6), effect = "logFC", significance = "adj p", contrast = "Condition")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("volcano accepts a non-syntactic significance column", {
  p <- .multigroup_volcano_labelled(
    .volcano_input(2),
    effect = "logFC",
    significance = "adj p",
    contrast = "Condition",
    label = "Name",
    yintercept = 0.5
  )
  built <- ggplot2::ggplot_build(p)
  expect_equal(sort(built$data[[1]]$y), sort(-log10(.volcano_input(2)[["adj p"]])))
})

test_that("plot_stat_violin_median() shows median and quartiles as one point range", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 20, with_missing = FALSE, seed = 1)
  lfq <- LFQData$new(istar$data, istar$config)
  stats <- lfq$get_Stats()$stats()
  p <- .plot_stat_violin_median(stats, "group_", stat = "meanAbundance")
  expect_length(p$layers, 2)
  summary_layer <- ggplot2::layer_data(p, 2)
  expect_false(anyNA(summary_layer$ymin))
  expect_true(all(summary_layer$ymin <= summary_layer$y & summary_layer$y <= summary_layer$ymax))
})

test_that("beeswarm point sizes stay within 1 to 4 mm for large child counts", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 5, with_missing = FALSE, seed = 1)
  lfq <- LFQData$new(istar$data, istar$config)
  pdata <- lfq$data_long()[lfq$data_long()$protein_Id == lfq$data_long()$protein_Id[1], ]
  pdata[[lfq$nr_children_col()]] <- rep(c(1, 40), length.out = nrow(pdata))
  p <- .plot_hierarchies_boxplot(pdata, "protein", lfq, beeswarm = TRUE, show_mean = FALSE)
  sizes <- ggplot2::layer_data(p, 2)$size
  expect_true(all(sizes >= 1 & sizes <= 4))
})
