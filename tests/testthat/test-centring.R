# Effect of centring on contrasts (srr RE2.3).

centring_fixture <- function() {
  sim <- sim_lfq_data_protein_config(n_proteins = 10, with_missing = FALSE, seed = 42)
  suppressMessages(LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq)
}

contrast_table <- function(lfq) {
  fit <- build_contrast_analysis(lfq, "~ group_", c(AvsC = "group_A - group_Ctrl", BvsC = "group_B - group_Ctrl"), "lm")
  dplyr::arrange(fit$get_contrasts(), .data$protein_Id, .data$contrast)
}

test_that("preserve_mean changes the average abundance, not the contrasts", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- centring_fixture()
  kept <- contrast_table(lfq$get_Transformer()$robscale(preserve_mean = TRUE)$lfq)
  zero <- contrast_table(lfq$get_Transformer()$robscale(preserve_mean = FALSE)$lfq)
  for (column in c("diff", "std.error", "statistic", "p.value")) {
    expect_equal(kept[[column]], zero[[column]], tolerance = 1e-10, label = column)
  }
  expect_gt(min(abs(kept$avgAbd - zero$avgAbd)), 1)
})

test_that("a shift common to all samples leaves contrasts unchanged; a per-group shift moves them by its difference", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- centring_fixture()
  base <- contrast_table(lfq)
  samples <- lfq$factors()
  common <- stats::setNames(rep(3, nrow(samples)), samples$sampleName)
  shifted <- contrast_table(lfq$get_Transformer()$subtract_sample_offsets(common)$lfq)
  expect_equal(shifted$diff, base$diff, tolerance = 1e-10)
  expect_equal(shifted$p.value, base$p.value, tolerance = 1e-10)
  expect_equal(shifted$avgAbd, base$avgAbd - 3, tolerance = 1e-10)
  by_group <- stats::setNames(ifelse(samples$group_ == "A", 0.5, 0), samples$sampleName)
  moved <- contrast_table(lfq$get_Transformer()$subtract_sample_offsets(by_group)$lfq)
  is_avsc <- moved$contrast == "AvsC"
  expect_equal(moved$diff[is_avsc], base$diff[is_avsc] - 0.5, tolerance = 1e-10)
  expect_equal(moved$diff[!is_avsc], base$diff[!is_avsc], tolerance = 1e-10)
})
