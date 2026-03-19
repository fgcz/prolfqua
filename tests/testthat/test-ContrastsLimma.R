test_that("strategy_limma returns correct structure", {
  strat <- prolfqua::strategy_limma("abundance ~ group_")
  expect_true(inherits(strat$formula, "formula"))
  expect_equal(strat$model_name, "limma")
  expect_false(strat$trend)
  expect_false(strat$robust)
})

test_that("build_model_limma returns ModelLimma", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)

  expect_true(inherits(mod, "ModelLimma"))
  expect_true(inherits(mod, "ModelInterface"))
  expect_true(!is.null(mod$fit))
  expect_true(!is.null(mod$dummy_model))
  expect_true(!is.null(mod$rowdata))
})

test_that("ModelLimma$get_coefficients works", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)

  coeffs <- mod$get_coefficients()
  expect_true(is.data.frame(coeffs))
  expect_true("Estimate" %in% colnames(coeffs))
  expect_true("Pr...t.." %in% colnames(coeffs))
  expect_true("factor" %in% colnames(coeffs))
  expect_true(nrow(coeffs) > 0)
})

test_that("ModelLimma$get_anova works", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)

  anova_tbl <- mod$get_anova()
  expect_true(is.data.frame(anova_tbl))
  expect_true("p.value" %in% colnames(anova_tbl))
  expect_true("F.value" %in% colnames(anova_tbl))
  expect_true(nrow(anova_tbl) > 0)
})

test_that("ContrastsLimma produces valid results", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)

  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contr_limma <- prolfqua::ContrastsLimma$new(mod, Contr)

  res <- contr_limma$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(all(
    c(
      "diff",
      "FDR",
      "p.value",
      "statistic",
      "std.error",
      "sigma",
      "df",
      "conf.low",
      "conf.high",
      "avgAbd",
      "modelName",
      "contrast"
    ) %in%
      colnames(res)
  ))
  expect_true(all(res$p.value >= 0 & res$p.value <= 1, na.rm = TRUE))
  expect_true(all(res$FDR >= 0 & res$FDR <= 1, na.rm = TRUE))
  # avgAbd can be NA for proteins with incomplete data (NA coefficients)
  expect_true(sum(!is.na(res$avgAbd)) > 0)
})

test_that("ContrastsLimma fold changes match prolfqua lm", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  # Limma pipeline
  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod_limma <- prolfqua::build_model_limma(lProt, strat)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contr_limma <- prolfqua::ContrastsLimma$new(mod_limma, Contr)
  res_limma <- contr_limma$get_contrasts()

  # prolfqua lm pipeline
  modelFunction <- prolfqua::strategy_lm("transformedIntensity ~ group_")
  mod_lm <- prolfqua::build_model(lProt, modelFunction)
  contr_lm <- prolfqua::Contrasts$new(mod_lm, Contr)
  res_lm <- contr_lm$get_contrasts()

  merged <- dplyr::inner_join(
    dplyr::select(res_limma, protein_Id, diff_limma = diff),
    dplyr::select(res_lm, protein_Id, diff_prolfqua = diff),
    by = "protein_Id"
  )

  # Fold changes should be essentially identical
  cc <- cor(merged$diff_limma, merged$diff_prolfqua, use = "complete.obs")
  expect_gt(cc, 0.99)
})

test_that("ContrastsLimma get_Plotter works", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contr_limma <- prolfqua::ContrastsLimma$new(mod, Contr)

  pl <- contr_limma$get_Plotter()
  expect_true(inherits(pl, "ContrastsPlotter"))
})

test_that("ContrastsLimma to_wide works", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contr_limma <- prolfqua::ContrastsLimma$new(mod, Contr)

  wide <- contr_limma$to_wide()
  expect_true(is.data.frame(wide))
  expect_true(ncol(wide) > 2)
})

test_that("ContrastsLimma get_contrast_sides works", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl", "B_vs_Ctrl" = "group_B - group_Ctrl")
  contr_limma <- prolfqua::ContrastsLimma$new(mod, Contr)

  sides <- contr_limma$get_contrast_sides()
  expect_equal(nrow(sides), 2)
  expect_true(all(c("group_1", "group_2") %in% colnames(sides)))
})

test_that("merge_contrasts_results works with ContrastsLimma", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.3)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contr_limma <- prolfqua::ContrastsLimma$new(mod, Contr)

  csi <- prolfqua::ContrastsMissing$new(lProt, contrasts = Contr)

  merged <- prolfqua::merge_contrasts_results(contr_limma, csi)
  expect_true(inherits(merged$merged, "ContrastsTable"))
  expect_true(nrow(merged$merged$get_contrasts()) > 0)
})

test_that("ContrastsLimma works with multiple contrasts", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)

  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl", "B_vs_Ctrl" = "group_B - group_Ctrl", "A_vs_B" = "group_A - group_B")
  contr_limma <- prolfqua::ContrastsLimma$new(mod, Contr)
  res <- contr_limma$get_contrasts()

  expect_equal(length(unique(res$contrast)), 3)
  # Each contrast should have the same number of proteins
  n_per_contrast <- table(res$contrast)
  expect_true(length(unique(n_per_contrast)) == 1)
})

test_that("build_model_limma passes weights from column name", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  # Add a weight column to the underlying data (per-sample weights)
  set.seed(42)
  sample_ids <- unique(lProt$data$sampleName)
  wt_map <- stats::setNames(runif(length(sample_ids), 0.1, 1.0), sample_ids)
  lProt$data$sample_weight <- wt_map[lProt$data$sampleName]

  # Fit with weights
  strat_wt <- prolfqua::strategy_limma("transformedIntensity ~ group_", weights = "sample_weight")
  mod_wt <- prolfqua::build_model_limma(lProt, strat_wt)

  # Fit without weights
  strat_nw <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod_nw <- prolfqua::build_model_limma(lProt, strat_nw)

  # Both should produce valid ModelLimma objects

  expect_true(inherits(mod_wt, "ModelLimma"))
  expect_true(inherits(mod_nw, "ModelLimma"))

  # Contrasts should work with weighted model
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contr_wt <- prolfqua::ContrastsLimma$new(mod_wt, Contr)
  res_wt <- contr_wt$get_contrasts()
  expect_true(is.data.frame(res_wt))
  expect_true(all(c("diff", "FDR", "p.value", "statistic") %in% colnames(res_wt)))

  contr_nw <- prolfqua::ContrastsLimma$new(mod_nw, Contr)
  res_nw <- contr_nw$get_contrasts()

  # Results should differ (weights are non-uniform)
  merged <- dplyr::inner_join(
    dplyr::select(res_wt, protein_Id, p_wt = p.value),
    dplyr::select(res_nw, protein_Id, p_nw = p.value),
    by = "protein_Id"
  )
  expect_false(all(abs(merged$p_wt - merged$p_nw) < 1e-10))
})

test_that("ContrastsLimma works with 2-factor design", {
  dd <- prolfqua::sim_lfq_data_2Factor_config(Nprot = 30, with_missing = FALSE)
  lProt <- prolfqua::LFQData$new(dd$data, dd$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ Treatment + Background")
  mod <- prolfqua::build_model_limma(lProt, strat)

  Contr <- c("A_vs_B" = "TreatmentA - TreatmentB", "X_vs_Z" = "BackgroundX - BackgroundZ")
  contr_limma <- prolfqua::ContrastsLimma$new(mod, Contr)
  res <- contr_limma$get_contrasts()

  expect_equal(length(unique(res$contrast)), 2)
  expect_true(all(c("diff", "FDR", "p.value") %in% colnames(res)))
})
