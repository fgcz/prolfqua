test_that("strategy_limma returns StrategyLimma R6 object", {
  strat <- prolfqua::strategy_limma("abundance ~ group_")
  expect_true(inherits(strat, "StrategyLimma"))
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

test_that("limma_impute contrasts use annotation matched to matrix columns", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 8)
  data <- istar$data |>
    dplyr::mutate(
      abundance = dplyr::case_when(
        .data$group_ == "A" ~ 20,
        .data$group_ == "Ctrl" ~ 17,
        TRUE ~ 12
      ) +
        as.numeric(factor(.data$protein_Id)) / 10 +
        as.numeric(factor(.data$sampleName)) / 100
    )
  sample_order <- unique(data$sampleName)
  interleaved <- sample_order[c(9, 2, 10, 1, 11, 3, 5, 12, 4, 6, 7, 8)]
  data <- data |>
    dplyr::mutate(.sample_order = match(.data$sampleName, interleaved)) |>
    dplyr::arrange(.data$.sample_order) |>
    dplyr::select(-".sample_order")

  lProt <- prolfqua::LFQData$new(data, istar$config)
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::ContrastsLimmaImputeFacade$new(lProt, "~ group_", contrasts, weights = NULL)
  res <- fa$get_contrasts()
  observed <- data |>
    dplyr::filter(.data$group_ %in% c("A", "Ctrl")) |>
    dplyr::summarise(mean_abundance = mean(.data$abundance), .by = c("protein_Id", "group_")) |>
    tidyr::pivot_wider(names_from = "group_", values_from = "mean_abundance") |>
    dplyr::mutate(expected_diff = .data$A - .data$Ctrl)

  comparison <- dplyr::inner_join(
    dplyr::select(res, "protein_Id", diff_limma = "diff"),
    dplyr::select(observed, "protein_Id", "expected_diff"),
    by = "protein_Id"
  )
  expect_equal(comparison$diff_limma, comparison$expected_diff, tolerance = 1e-8)
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
  sample_ids <- unique(lProt$data_long()$sampleName)
  wt_map <- stats::setNames(runif(length(sample_ids), 0.1, 1.0), sample_ids)
  lProt$set_data(
    lProt$data_long() |> dplyr::mutate(sample_weight = wt_map[sampleName])
  )

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
  dd <- prolfqua::sim_lfq_data_2factor_config(Nprot = 30, with_missing = FALSE)
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


# --- build_model_limma_impute tests ---

test_that("build_model_limma_impute recovers failed proteins", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")

  mod_na <- prolfqua::build_model_limma(lfqdata, strat)
  mod_imp <- prolfqua::build_model_limma_impute(lfqdata, strat)

  na_count_plain <- sum(rowSums(is.na(mod_na$fit$coefficients)) > 0)
  na_count_imputed <- sum(rowSums(is.na(mod_imp$fit$coefficients)) > 0)

  # Imputed model should have fewer (ideally zero) NA coefficient rows

  expect_true(na_count_imputed < na_count_plain || na_count_plain == 0)
  expect_equal(na_count_imputed, 0)
})

test_that("ContrastsLimmaImputeFacade end-to-end", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::ContrastsLimmaImputeFacade$new(lfqdata, "~ group_", contrasts)
  res <- fa$get_contrasts()
  expect_true(nrow(res) > 0)
  expect_true(all(c("diff", "FDR", "p.value", "statistic") %in% colnames(res)))

  # get_Plotter returns a ContrastsPlotter
  pl <- fa$get_Plotter()
  expect_true(inherits(pl, "ContrastsPlotter"))

  # to_wide works
  wide <- fa$to_wide()
  expect_true(nrow(wide) > 0)
})

test_that("limma_impute fold changes match plain limma for non-imputed proteins", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.3)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa_plain <- prolfqua::ContrastsLimmaFacade$new(lfqdata, "~ group_", contrasts)
  fa_imp <- prolfqua::ContrastsLimmaImputeFacade$new(lfqdata, "~ group_", contrasts)

  res_plain <- fa_plain$get_contrasts()
  res_imp <- fa_imp$get_contrasts()

  # Imputed should have at least as many rows
  expect_true(nrow(res_imp) >= nrow(res_plain))

  # Fold changes for shared proteins should correlate highly
  merged <- dplyr::inner_join(
    dplyr::select(res_plain, protein_Id, diff_plain = diff),
    dplyr::select(res_imp, protein_Id, diff_imp = diff),
    by = "protein_Id"
  )
  if (nrow(merged) >= 3) {
    expect_true(cor(merged$diff_plain, merged$diff_imp, use = "complete.obs") > 0.95)
  }
})

test_that("build_contrast_analysis dispatches limma_impute", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20, weight_missing = 0.3)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma_impute")
  expect_true(inherits(fa, "ContrastsLimmaImputeFacade"))
  res <- fa$get_contrasts()
  expect_true(nrow(res) > 0)
})

test_that("df_method borrowed works for limma_impute", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::ContrastsLimmaImputeFacade$new(
    lfqdata,
    "~ group_",
    contrasts,
    df_method = "borrowed"
  )
  res <- fa$get_contrasts()
  expect_true(nrow(res) > 0)
})

test_that("df correction: imputed proteins do not use inflated df", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")

  mod_na <- prolfqua::build_model_limma(lfqdata, strat)
  mod_imp <- prolfqua::build_model_limma_impute(lfqdata, strat)

  # Identify which proteins were imputed (had NA coefficients in plain fit)
  failed <- which(rowSums(is.na(mod_na$fit$coefficients)) > 0)

  if (length(failed) > 0) {
    p <- ncol(mod_imp$fit$coefficients)
    wide <- lfqdata$data_wide(as.matrix = TRUE)
    n_observed <- rowSums(!is.na(wide$data))
    max_df_full <- ncol(wide$data) - p

    # Imputed proteins should have df < max possible (not counting imputed as real)
    imputed_df <- mod_imp$fit$df.residual[failed]
    expect_true(all(imputed_df <= max_df_full))
    # df should be based on n_observed, not n_total
    expected_df <- pmax(n_observed[failed] - p, 1)
    expect_equal(unname(imputed_df), unname(expected_df))
  }
})


# ---- limma_voom tests ----

test_that("build_model_limma_voom returns ModelLimma", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma_voom(lProt, strat)

  expect_true(inherits(mod, "ModelLimma"))
  expect_true(inherits(mod, "ModelInterface"))
  expect_true(!is.null(mod$fit))
  expect_true(!is.null(mod$dummy_model))
})

test_that("build_model_limma_voom produces different weights than plain limma", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa_limma <- prolfqua::ContrastsLimmaFacade$new(lProt, "~ group_", contrasts)
  fa_voom <- prolfqua::ContrastsLimmaVoomFacade$new(lProt, "~ group_", contrasts)

  res_limma <- fa_limma$get_contrasts()
  res_voom <- fa_voom$get_contrasts()

  # Both should produce results

  expect_true(nrow(res_voom) > 0)
  expect_true(nrow(res_limma) > 0)

  # Fold changes should be correlated but not identical (weights differ)
  merged <- dplyr::inner_join(
    dplyr::select(res_limma, protein_Id, diff_limma = diff),
    dplyr::select(res_voom, protein_Id, diff_voom = diff),
    by = "protein_Id"
  )
  expect_true(cor(merged$diff_limma, merged$diff_voom, use = "complete.obs") > 0.9)
})

test_that("build_contrast_analysis dispatches limma_voom", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma_voom")
  expect_true(inherits(fa, "ContrastsLimmaVoomFacade"))
  res <- fa$get_contrasts()
  expect_true(nrow(res) > 0)
  expect_equal(unique(res$modelName), "limma_voom")
})

test_that("build_contrast_analysis dispatches limma_voom_impute", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma_voom_impute")
  expect_true(inherits(fa, "ContrastsLimmaVoomImputeFacade"))
  res <- fa$get_contrasts()
  expect_true(nrow(res) > 0)
  expect_equal(unique(res$modelName), "limma_voom_impute")
})

test_that("limma_voom_impute recovers more proteins than limma_voom", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa_voom <- prolfqua::ContrastsLimmaVoomFacade$new(lfqdata, "~ group_", contrasts)
  fa_voom_imp <- prolfqua::ContrastsLimmaVoomImputeFacade$new(lfqdata, "~ group_", contrasts)

  res_voom <- fa_voom$get_contrasts()
  res_imp <- fa_voom_imp$get_contrasts()

  expect_true(nrow(res_imp) >= nrow(res_voom))
})

test_that("build_contrast_analysis dispatches deqms_voom", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "deqms_voom")
  expect_true(inherits(fa, "ContrastsDEqMSVoomFacade"))
  res <- fa$get_contrasts()
  expect_true(nrow(res) > 0)
  expect_equal(unique(res$modelName), "deqms_voom")
})

test_that("limma_voom with external weights (nr_children)", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  # With default weights (nr_children)
  fa_wt <- prolfqua::ContrastsLimmaVoomFacade$new(lfqdata, "~ group_", contrasts)
  # Without weights
  fa_nowt <- prolfqua::ContrastsLimmaVoomFacade$new(lfqdata, "~ group_", contrasts, weights = NULL)

  res_wt <- fa_wt$get_contrasts()
  res_nowt <- fa_nowt$get_contrasts()

  expect_true(nrow(res_wt) > 0)
  expect_true(nrow(res_nowt) > 0)
})


test_that("limma_impute flags LOD-rescued proteins via estimate_type", {
  # Same fixture as the lm_impute regression test in test-ImputeModel.R;
  # weight_missing = 0.5 guarantees some proteins have NA limma
  # coefficients, which build_model_limma_impute then refits at LOD.
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5, seed = 42)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  # Plain limma: no imputed_proteins, modelName is the facade key, all observed.
  fa_plain <- prolfqua::ContrastsLimmaFacade$new(lfqdata, "~ group_", contrasts)
  res_plain <- fa_plain$get_contrasts()
  expect_setequal(unique(res_plain$modelName), "limma")
  expect_setequal(unique(res_plain$estimate_type), "observed")
  expect_length(fa_plain$model$imputed_proteins, 0)

  # limma_impute: modelName is the facade key; refit proteins are flagged
  # in estimate_type as "lod_imputed".
  fa_imp <- prolfqua::ContrastsLimmaImputeFacade$new(lfqdata, "~ group_", contrasts)
  expect_gt(length(fa_imp$model$imputed_proteins), 0)
  res_imp <- fa_imp$get_contrasts()
  expect_setequal(unique(res_imp$modelName), "limma_impute")
  expect_setequal(unique(res_imp$estimate_type), c("observed", "lod_imputed"))
  # Refit count should match imputed_proteins (one row per contrast,
  # and we asked for one contrast).
  expect_equal(
    sum(res_imp$estimate_type == "lod_imputed"),
    length(fa_imp$model$imputed_proteins)
  )
})


test_that("limma_voom_impute also flags LOD-rescued proteins", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 50, weight_missing = 0.5, seed = 42)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::ContrastsLimmaVoomImputeFacade$new(lfqdata, "~ group_", contrasts)
  expect_gt(length(fa$model$imputed_proteins), 0)
  res <- fa$get_contrasts()
  expect_setequal(unique(res$modelName), "limma_voom_impute")
  expect_true("lod_imputed" %in% res$estimate_type)
})

test_that("limma_voom_impute returns an ordinary fit when no rows fail", {
  istar <- prolfqua::sim_lfq_data_protein_config(
    Nprot = 30,
    weight_missing = 0,
    seed = 42
  )
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  strategy <- prolfqua::strategy_limma(
    "transformedIntensity ~ group_"
  )

  model <- prolfqua::build_model_limma_voom_impute(
    lfqdata,
    strategy
  )

  expect_s3_class(model, "ModelLimma")
  expect_length(model$imputed_proteins, 0)
  expect_false(anyNA(model$fit$coefficients))
})

test_that("limma_voom_impute supports borrowed residual degrees of freedom", {
  istar <- prolfqua::sim_lfq_data_protein_config(
    Nprot = 50,
    weight_missing = 0.5,
    seed = 42
  )
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  strategy <- prolfqua::strategy_limma(
    "transformedIntensity ~ group_"
  )

  model <- prolfqua::build_model_limma_voom_impute(
    lfqdata,
    strategy,
    df_method = "borrowed"
  )
  imputed_rows <- model$rowdata$protein_Id %in%
    model$imputed_proteins

  expect_true(any(imputed_rows))
  expect_length(
    unique(model$fit$df.residual[imputed_rows]),
    1
  )
})


test_that("ContrastsLMMissingFacade emits a deprecation warning", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  expect_warning(
    prolfqua::ContrastsLMMissingFacade$new(lfqdata, "~ group_", contrasts),
    "deprecated"
  )
})


test_that("ContrastsMissing emits a deprecation warning at construction", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  expect_warning(
    prolfqua::ContrastsMissing$new(lfqdata, contrasts = contrasts),
    "deprecated"
  )
})
