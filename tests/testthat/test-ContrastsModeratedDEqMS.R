# Test ContrastsModeratedDEqMS R6 class

test_that("find_d0_deqms returns sensible values", {
  # Positive input: should return finite d0
  d0 <- find_d0_deqms(0.5)
  expect_true(is.finite(d0))
  expect_true(d0 > 0)

  # Zero or negative: should return Inf
  expect_equal(find_d0_deqms(0), Inf)
  expect_equal(find_d0_deqms(-1), Inf)

  # NA: should return Inf
  expect_equal(find_d0_deqms(NA), Inf)

  # Very small positive: large d0 (strong prior)
  d0_small <- find_d0_deqms(0.001)
  d0_large <- find_d0_deqms(1.0)
  expect_true(d0_small > d0_large)
})

test_that("moderated_p_deqms works on simulated data", {
  istar <- sim_lfq_data_protein_config(Nprot = 50)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  modelFunction <- strategy_lm("transformedIntensity ~ group_")
  mod <- build_model(lProt, modelFunction)

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  contrast <- Contrasts$new(mod, Contr)
  contrast_result <- contrast$get_contrasts(all = FALSE)

  # Add count column
  count_df <- dplyr::select(istar$data, dplyr::all_of(c(istar$config$hierarchy_keys_depth(), "nr_peptides"))) |>
    dplyr::distinct()
  contrast_result <- dplyr::inner_join(contrast_result, count_df, by = istar$config$hierarchy_keys_depth())
  contrast_result$nr_peptides <- pmax(contrast_result$nr_peptides, 1)

  # Test single-group moderation
  result <- moderated_p_deqms(contrast_result, count_col = "nr_peptides")

  expect_true("moderated.statistic" %in% colnames(result))
  expect_true("moderated.p.value" %in% colnames(result))
  expect_true("moderated.var.post" %in% colnames(result))
  expect_true("moderated.var.prior" %in% colnames(result))
  expect_true("moderated.df.prior" %in% colnames(result))
  expect_true("moderated.df.total" %in% colnames(result))
  expect_true("moderated.conf.low" %in% colnames(result))
  expect_true("moderated.conf.high" %in% colnames(result))

  # All p-values should be in (0, 1]
  expect_true(all(result$moderated.p.value > 0))
  expect_true(all(result$moderated.p.value <= 1))

  # Prior variance should vary per protein (count-dependent)
  expect_true(length(unique(result$moderated.var.prior)) > 1)
})

test_that("ContrastsModeratedDEqMS works end-to-end", {
  istar <- sim_lfq_data_protein_config(Nprot = 50)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  modelFunction <- strategy_lm("transformedIntensity ~ group_")
  mod <- build_model(lProt, modelFunction)

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  contrast <- Contrasts$new(mod, Contr)

  # Build count_df from data
  count_df <- dplyr::select(istar$data, dplyr::all_of(c(istar$config$hierarchy_keys_depth(), "nr_peptides"))) |>
    dplyr::distinct()

  deqms <- ContrastsModeratedDEqMS$new(
    contrast,
    count_df = count_df,
    count_column = "nr_peptides"
  )

  # get_contrasts returns expected columns
  x <- deqms$get_contrasts()
  expect_s3_class(x, "data.frame")
  expect_true(all(c("diff", "p.value", "FDR", "sigma", "statistic", "df", "conf.low", "conf.high") %in% colnames(x)))
  expect_true(all(x$p.value > 0 & x$p.value <= 1))
  expect_true(all(x$FDR >= 0 & x$FDR <= 1))

  # modelName should end with _DEqMS
  expect_true(all(grepl("_DEqMS", x$modelName)))

  # get_contrasts(all = TRUE) returns both original and moderated columns
  x_all <- deqms$get_contrasts(all = TRUE)
  expect_true("moderated.statistic" %in% colnames(x_all))
  expect_true("moderated.p.value" %in% colnames(x_all))
  expect_true("FDR.moderated" %in% colnames(x_all))

  # get_contrast_sides delegates correctly
  cs <- deqms$get_contrast_sides()
  expect_s3_class(cs, "data.frame")
  expect_true(all(c("group_1", "group_2") %in% colnames(cs)))

  # get_linfct delegates
  lf <- deqms$get_linfct()
  expect_true(!is.null(lf))

  # to_wide returns wide format
  tw <- deqms$to_wide()
  expect_s3_class(tw, "data.frame")

  # get_Plotter returns ContrastsPlotter
  pl <- deqms$get_Plotter()
  expect_s3_class(pl, "ContrastsPlotter")
})

test_that("ContrastsModeratedDEqMS works with multiple contrasts", {
  istar <- sim_lfq_data_protein_config(Nprot = 30)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  modelFunction <- strategy_lm("transformedIntensity ~ group_")
  mod <- build_model(lProt, modelFunction)

  Contr <- c(
    "AvsCtrl" = "group_A - group_Ctrl",
    "BvsCtrl" = "group_B - group_Ctrl"
  )
  contrast <- Contrasts$new(mod, Contr)

  count_df <- dplyr::select(istar$data, dplyr::all_of(c(istar$config$hierarchy_keys_depth(), "nr_peptides"))) |>
    dplyr::distinct()

  deqms <- ContrastsModeratedDEqMS$new(
    contrast,
    count_df = count_df,
    count_column = "nr_peptides"
  )

  x <- deqms$get_contrasts()
  expect_true(all(c("AvsCtrl", "BvsCtrl") %in% x$contrast))
})

test_that("merge_contrasts_results works with ContrastsModeratedDEqMS", {
  istar <- sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.4)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")

  modelFunction <- strategy_lm("transformedIntensity ~ group_")
  mod <- build_model(lProt, modelFunction)
  contrast <- Contrasts$new(mod, Contr)

  count_df <- dplyr::select(istar$data, dplyr::all_of(c(istar$config$hierarchy_keys_depth(), "nr_peptides"))) |>
    dplyr::distinct()

  deqms <- ContrastsModeratedDEqMS$new(
    contrast,
    count_df = count_df,
    count_column = "nr_peptides"
  )

  csi <- ContrastsMissing$new(lProt, contrasts = Contr)

  merged <- merge_contrasts_results(deqms, csi)
  expect_true("merged" %in% names(merged))
  expect_true("more" %in% names(merged))
  expect_true("same" %in% names(merged))

  merged_contrasts <- merged$merged$get_contrasts()
  expect_s3_class(merged_contrasts, "data.frame")
})

test_that("ContrastsModeratedDEqMS volcano renders", {
  istar <- sim_lfq_data_protein_config(Nprot = 30)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  modelFunction <- strategy_lm("transformedIntensity ~ group_")
  mod <- build_model(lProt, modelFunction)

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  contrast <- Contrasts$new(mod, Contr)

  count_df <- dplyr::select(istar$data, dplyr::all_of(c(istar$config$hierarchy_keys_depth(), "nr_peptides"))) |>
    dplyr::distinct()

  deqms <- ContrastsModeratedDEqMS$new(
    contrast,
    count_df = count_df,
    count_column = "nr_peptides"
  )

  pl <- deqms$get_Plotter()
  vp <- pl$volcano()
  expect_true(!is.null(vp))
})
