test_that("strategy_limma returns StrategyLimma R6 object", {
  strat <- prolfqua::strategy_limma("abundance ~ group_")
  expect_true(inherits(strat, "StrategyLimma"))
  expect_true(inherits(strat$formula, "formula"))
  expect_equal(strat$model_name, "limma")
  expect_false(strat$trend)
  expect_false(strat$robust)
})

test_that("build_model_limma returns ModelLimma", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)

  expect_true(inherits(mod, "ModelLimma"))
  expect_true(inherits(mod, "ModelInterface"))
  # one fit row per protein, one coefficient per level of group_
  expect_s4_class(mod$fit, "MArrayLM")
  expect_equal(dim(mod$fit$coefficients), c(20L, 3L))
  expect_equal(names(stats::coef(mod$dummy_model)), colnames(mod$fit$coefficients))
  expect_equal(mod$rowdata$protein_Id, rownames(mod$fit$coefficients))
})

test_that("ModelLimma$get_coefficients works", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)

  coeffs <- mod$get_coefficients()
  expect_true(is.data.frame(coeffs))
  expect_true("Estimate" %in% colnames(coeffs))
  expect_true("Pr...t.." %in% colnames(coeffs))
  expect_true("factor" %in% colnames(coeffs))
  # 20 proteins x 3 coefficients
  expect_equal(nrow(coeffs), 20 * 3)
  expect_setequal(coeffs$factor, c("(Intercept)", "group_B", "group_Ctrl"))
})

test_that("ModelLimma$get_anova works", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)

  anova_tbl <- mod$get_anova()
  expect_true(is.data.frame(anova_tbl))
  expect_true("p.value" %in% colnames(anova_tbl))
  expect_true("F.value" %in% colnames(anova_tbl))
  # one F-test of the group_ coefficients per protein
  expect_equal(nrow(anova_tbl), 20)
  expect_setequal(anova_tbl$factor, "group_B+group_Ctrl")
})

test_that("ContrastsLimma produces valid results", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0)
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
      "std.error.unmoderated",
      "sigma",
      "df",
      "df.unmoderated",
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
  expect_equal(
    unname(res$statistic),
    unname(res$diff / res$std.error),
    tolerance = 1e-12
  )
  expect_true(any(abs(res$std.error - res$std.error.unmoderated) > 1e-12, na.rm = TRUE))
  expect_true(any(abs(res$df - res$df.unmoderated) > 1e-12, na.rm = TRUE))
})

test_that("ContrastsLimma exposes identical standard and raw pairs without eBayes", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")
  mod <- prolfqua::build_model_limma(
    lProt,
    prolfqua::strategy_limma("transformedIntensity ~ group_")
  )
  contrast <- prolfqua::ContrastsLimma$new(
    mod,
    c("A_vs_Ctrl" = "group_A - group_Ctrl"),
    eBayes = FALSE
  )

  res <- contrast$get_contrasts()

  expect_equal(res$std.error.unmoderated, res$std.error)
  expect_equal(res$df.unmoderated, res$df)
})

test_that("ContrastsLimma fold changes match prolfqua lm", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 8)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 20)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.3)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma(lProt, strat)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contr_limma <- prolfqua::ContrastsLimma$new(mod, Contr)

  imp <- prolfqua::ContrastsLMImputeFacade$new(lProt, "~ group_", Contr)$contrast

  merged <- prolfqua::merge_contrasts_results(contr_limma, imp)
  expect_true(inherits(merged$merged, "ContrastsTable"))
  # limma cannot estimate a protein without group A observations; the imputed lm analysis adds it
  merged_contrasts <- merged$merged$get_contrasts()
  expect_equal(nrow(merged_contrasts), dplyr::n_distinct(lProt$data_long()$protein_Id))
  expect_equal(nrow(merged$more$get_contrasts()), sum(is.na(contr_limma$get_contrasts()$statistic)))
})

test_that("ContrastsLimma works with multiple contrasts", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  # Add a weight column to the underlying data (per-sample weights)
  withr::local_seed(42)
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
  dd <- prolfqua::sim_lfq_data_2factor_config(n_proteins = 30, with_missing = FALSE)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.5)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0.5)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::ContrastsLimmaImputeFacade$new(lfqdata, "~ group_", contrasts)
  res <- fa$get_contrasts()
  # imputation gives every protein an estimate
  n_proteins <- dplyr::n_distinct(lfqdata$data_long()$protein_Id)
  expect_equal(nrow(res), n_proteins)
  expect_true(all(c("diff", "FDR", "p.value", "statistic") %in% colnames(res)))

  # get_Plotter returns a ContrastsPlotter
  pl <- fa$get_Plotter()
  expect_true(inherits(pl, "ContrastsPlotter"))

  # to_wide works
  wide <- fa$to_wide()
  expect_equal(nrow(wide), n_proteins)
})

test_that("limma_impute fold changes match plain limma for non-imputed proteins", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.3)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 20, weight_missing = 0.3)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma_impute")
  expect_true(inherits(fa, "ContrastsLimmaImputeFacade"))
  res <- fa$get_contrasts()
  expect_equal(nrow(res), dplyr::n_distinct(lfqdata$data_long()$protein_Id))
})

test_that("df_method borrowed works for limma_impute", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0.5)
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
  expect_equal(nrow(res), dplyr::n_distinct(lfqdata$data_long()$protein_Id))
})

test_that("df correction: imputed proteins do not use inflated df", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.5)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 20)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  strat <- prolfqua::strategy_limma("transformedIntensity ~ group_")
  mod <- prolfqua::build_model_limma_voom(lProt, strat)

  expect_true(inherits(mod, "ModelLimma"))
  expect_true(inherits(mod, "ModelInterface"))
  # one fit row and one row of vooma weights per protein
  expect_s4_class(mod$fit, "MArrayLM")
  expect_equal(dim(mod$fit$coefficients), c(20L, 3L))
  expect_equal(dim(mod$weights), c(20L, dplyr::n_distinct(lProt$data_long()[[lProt$sample_name()]])))
  expect_equal(names(stats::coef(mod$dummy_model)), colnames(mod$fit$coefficients))
})

# complete protein data and the limma::vooma() pipeline that the voom backends must reproduce
vooma_reference <- function() {
  sim <- prolfqua::sim_lfq_data_protein_config(
    n_proteins = 200,
    weight_missing = 0,
    with_missing = FALSE,
    seed = 5
  )
  lfqdata <- prolfqua::LFQData$new(sim$data, sim$config)
  lfqdata$rename_response("transformedIntensity")
  wide <- lfqdata$data_wide(as.matrix = TRUE)
  design <- stats::model.matrix(~group_, data = wide$annotation)
  voom <- limma::vooma(wide$data, design, span = 0.5)
  fit <- limma::lmFit(voom, design)
  # A - Ctrl with A as the intercept level
  contrast_matrix <- matrix(c(0, 0, -1), ncol = 1, dimnames = list(colnames(design), "A_vs_Ctrl"))
  top <- limma::topTable(
    limma::eBayes(limma::contrasts.fit(fit, contrast_matrix)),
    number = Inf,
    sort.by = "none"
  )
  list(sim = sim, lfqdata = lfqdata, design = design, voom = voom, fit = fit, top = top)
}

test_that("limma_voom and limma_voom_impute equal limma::vooma() on complete data", {
  ref <- vooma_reference()
  expect_identical(colnames(ref$design), c("(Intercept)", "group_B", "group_Ctrl"))
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  for (facade in list(prolfqua::ContrastsLimmaVoomFacade, prolfqua::ContrastsLimmaVoomImputeFacade)) {
    fa <- facade$new(ref$lfqdata, "~ group_", contrasts, weights = NULL)
    expect_equal(unname(fa$model$weights), unname(ref$voom$weights))
    expect_equal(unname(fa$model$fit$stdev.unscaled), unname(ref$fit$stdev.unscaled))
    res <- fa$get_contrasts()
    expect_equal(res$diff, ref$top[res$protein_Id, "logFC"])
    expect_equal(res$p.value, ref$top[res$protein_Id, "P.Value"])
  }
})

test_that("deqms_voom fits with the limma::vooma() weights", {
  ref <- vooma_reference()
  fa <- prolfqua::build_contrast_analysis(
    ref$lfqdata,
    "~ group_",
    c("A_vs_Ctrl" = "group_A - group_Ctrl"),
    method = "deqms_voom"
  )
  expect_equal(unname(fa$model$weights), unname(ref$voom$weights))
  expect_equal(unname(fa$model$fit$stdev.unscaled), unname(ref$fit$stdev.unscaled))
})

test_that("vooma weights change p-values but not one-factor fold changes", {
  ref <- vooma_reference()
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  res_limma <- prolfqua::ContrastsLimmaFacade$new(ref$lfqdata, "~ group_", contrasts, weights = NULL)$get_contrasts()
  res_voom <- prolfqua::ContrastsLimmaVoomFacade$new(ref$lfqdata, "~ group_", contrasts, weights = NULL)$get_contrasts()
  res_voom <- res_voom[match(res_limma$protein_Id, res_voom$protein_Id), ]

  # all samples of a group share one weight, so the group means stay unweighted means
  expect_equal(res_voom$diff, res_limma$diff)
  expect_gt(max(abs(res_voom$p.value - res_limma$p.value)), 1e-3)
})

test_that("build_contrast_analysis dispatches limma_voom", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma_voom")
  expect_true(inherits(fa, "ContrastsLimmaVoomFacade"))
  res <- fa$get_contrasts()
  # every protein is observed in groups A and Ctrl
  expect_equal(nrow(res), dplyr::n_distinct(lfqdata$data_long()$protein_Id))
  expect_equal(unique(res$modelName), "limma_voom")
})

test_that("build_contrast_analysis dispatches limma_voom_impute", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0.5)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limma_voom_impute")
  expect_true(inherits(fa, "ContrastsLimmaVoomImputeFacade"))
  res <- fa$get_contrasts()
  expect_equal(nrow(res), dplyr::n_distinct(lfqdata$data_long()$protein_Id))
  expect_equal(unique(res$modelName), "limma_voom_impute")
})

test_that("limma_voom_impute recovers more proteins than limma_voom", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.5)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "deqms_voom")
  expect_true(inherits(fa, "ContrastsDEqMSVoomFacade"))
  expect_no_warning(res <- fa$get_contrasts())
  # one row per protein that can estimate the contrast; get_missing() lists the others (an empty group)
  expect_equal(nrow(res) + nrow(fa$get_missing()), 30)
  expect_false(anyNA(res$diff))
  expect_equal(unique(res$modelName), "deqms_voom")
  # proteins without residual degrees of freedom get NA p-values, all others finite ones
  no_df <- res$df.unmoderated == 0
  expect_true(all(is.na(res$p.value[no_df])))
  expect_true(all(is.finite(res$p.value[!no_df])))
})

test_that("limma_voom external weights (nr_children) act as relative precisions", {
  ref <- vooma_reference()
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  res_none <- prolfqua::ContrastsLimmaVoomFacade$new(ref$lfqdata, "~ group_", contrasts, weights = NULL)$get_contrasts()

  # a constant weight rescales the preliminary residual variance and the trend by the same factor
  constant <- ref$sim
  constant$data[[ref$lfqdata$nr_children_col()]] <- 3
  lfq_constant <- prolfqua::LFQData$new(constant$data, constant$config)
  lfq_constant$rename_response("transformedIntensity")
  res_constant <- prolfqua::ContrastsLimmaVoomFacade$new(lfq_constant, "~ group_", contrasts)$get_contrasts()
  res_constant <- res_constant[match(res_none$protein_Id, res_constant$protein_Id), ]
  expect_equal(res_constant$p.value, res_none$p.value)

  # varying peptide counts per protein and sample change the fit
  res_counts <- prolfqua::ContrastsLimmaVoomFacade$new(ref$lfqdata, "~ group_", contrasts)$get_contrasts()
  res_counts <- res_counts[match(res_none$protein_Id, res_counts$protein_Id), ]
  expect_gt(max(abs(res_counts$p.value - res_none$p.value)), 1e-3)
})


test_that("limma_impute flags LOD-rescued proteins via estimate_type", {
  # Same fixture as the lm_impute regression test in test-ImputeModel.R;
  # weight_missing = 0.5 guarantees some proteins have NA limma
  # coefficients, which build_model_limma_impute then refits at LOD.
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.5, seed = 42)
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
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.5, seed = 42)
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
    n_proteins = 30,
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
    n_proteins = 50,
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

test_that("matrix models retain the inputs and weights used for their coefficients", {
  sim <- sim_lfq_data_protein_config(n_proteins = 20, weight_missing = 0, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)
  for (method in c("limma", "limma_voom")) {
    facade <- build_contrast_analysis(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"), method)
    model <- facade$model
    refit <- limma::lmFit(model$E, model$fit$design, weights = model$weights)
    expect_equal(refit$coefficients, model$fit$coefficients)
    expect_equal(unname(rowSums(model$observed)), rep(ncol(model$E), nrow(model$E)))
    expect_equal(model$annotation[[lfq$sample_name()]], colnames(model$E))
    expect_equal(model$formula, stats::as.formula(paste(lfq$response(), "~ group_")))
  }
})

test_that("limma imputation retains rescued responses and the original observation mask", {
  sim <- sim_lfq_data_protein_config(n_proteins = 20, weight_missing = 0, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)
  dat <- lfq$data_long()
  protein <- unique(dat$protein_Id)[1]
  lfq$set_data(dat[!(dat$protein_Id == protein & dat$group_ == "A"), ])
  original <- lfq$data_wide(as.matrix = TRUE)$data
  for (method in c("limma_impute", "limma_voom_impute")) {
    facade <- build_contrast_analysis(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"), method)
    model <- facade$model
    row <- which(model$rowdata$protein_Id == protein)
    expect_equal(model$observed, is.finite(original))
    expect_equal(model$E[-row, ], original[-row, ])
    expect_true(all(is.finite(model$E[row, ])))
    refit <- limma::lmFit(model$E, model$fit$design, weights = model$weights)
    expect_equal(refit$coefficients, model$fit$coefficients)
  }
})

test_that("facades retain a snapshot independent of the caller", {
  sim <- sim_lfq_data_protein_config(n_proteins = 10, weight_missing = 0, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)
  facade <- build_contrast_analysis(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"))
  before <- facade$.lfqdata$data_long()
  dat <- lfq$data_long()
  dat[[lfq$response()]] <- dat[[lfq$response()]] + 100
  lfq$set_data(dat)
  config <- lfq$get_config()
  config$sample_name <- "changed"
  expect_equal(facade$.lfqdata$data_long(), before)
  expect_false(identical(facade$.lfqdata$sample_name(), "changed"))
})

test_that("build_model_limma weights from a matrix or an annotation column equal limma::lmFit", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0, seed = 5)
  lProt <- prolfqua::LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")
  wide <- lProt$data_wide(as.matrix = TRUE)
  design <- model.matrix(~group_, data = wide$annotation)

  # observation weights as a matrix
  W <- matrix(
    withr::with_seed(1, stats::runif(length(wide$data), 0.2, 1)),
    nrow = nrow(wide$data),
    dimnames = dimnames(wide$data)
  )
  mod <- build_model_limma(lProt, strategy_limma("transformedIntensity ~ group_", weights = W))
  ref <- limma::lmFit(wide$data, design, weights = W)
  expect_equal(unname(mod$fit$coefficients), unname(ref$coefficients))
  expect_equal(unname(mod$fit$stdev.unscaled), unname(ref$stdev.unscaled))

  # per-sample weights from a column of the sample annotation
  w_sample <- stats::setNames(
    withr::with_seed(2, stats::runif(nrow(wide$annotation), 0.2, 1)),
    wide$annotation$sampleName
  )
  data <- lProt$data_long()
  data$sample_weight <- w_sample[data$sampleName]
  config <- lProt$get_config()$clone(deep = TRUE)
  config$factors[["sample_weight"]] <- "sample_weight"
  lw <- LFQData$new(data, config)
  wide_w <- lw$data_wide(as.matrix = TRUE)
  expect_true("sample_weight" %in% colnames(wide_w$annotation))
  mod_w <- build_model_limma(lw, strategy_limma("transformedIntensity ~ group_", weights = "sample_weight"))
  ref_w <- limma::lmFit(
    wide_w$data,
    model.matrix(~group_, data = wide_w$annotation),
    weights = wide_w$annotation$sample_weight
  )
  expect_equal(unname(mod_w$fit$coefficients), unname(ref_w$coefficients))
  expect_equal(unname(mod_w$fit$stdev.unscaled), unname(ref_w$stdev.unscaled))
})

# --- proteins with an empty group or cell: only the estimable contrasts, equal to the observed means ---

# Differences and get_missing() contrasts of one protein, for one method
protein_contrasts <- function(lfq, modelstr, contrasts, method, protein) {
  fa <- suppressMessages(prolfqua::build_contrast_analysis(lfq, modelstr, contrasts, method = method))
  # deqms_voom's variance-by-count loess warns on these small simulations
  res <- suppressWarnings(suppressMessages(fa$get_contrasts()))
  res <- res[res$protein_Id == protein, ]
  missing <- suppressWarnings(suppressMessages(fa$get_missing()))
  list(
    diff = stats::setNames(res$diff, res$contrast),
    missing = sort(missing$contrast[missing$protein_Id == protein])
  )
}

one_factor_contrasts <- c(
  A_vs_Ctrl = "group_A - group_Ctrl",
  B_vs_Ctrl = "group_B - group_Ctrl",
  B_vs_A = "group_B - group_A"
)

test_that("a protein without its reference group gets only its estimable contrasts in every method", {
  sim <- prolfqua::sim_lfq_data_protein_config(n_proteins = 40, weight_missing = 0, seed = 1)
  protein <- sort(unique(sim$data$protein_Id))[1]
  data <- sim$data[!(sim$data$protein_Id == protein & sim$data$group_ == "A"), ]
  lfq <- prolfqua::LFQData$new(data, sim$config)
  observed <- data[data$protein_Id == protein, ]
  means <- tapply(observed$abundance, observed$group_, mean)

  for (method in c("lm", "limma", "limma_voom", "deqms_voom")) {
    res <- protein_contrasts(lfq, "~ group_", one_factor_contrasts, method, protein)
    expect_equal(res$diff, c(B_vs_Ctrl = unname(means["B"] - means["Ctrl"])), info = method)
    expect_equal(res$missing, c("A_vs_Ctrl", "B_vs_A"), info = method)
  }
})

test_that("a protein without a non-reference group gets only its estimable contrasts in every method", {
  sim <- prolfqua::sim_lfq_data_protein_config(n_proteins = 40, weight_missing = 0, seed = 1)
  protein <- sort(unique(sim$data$protein_Id))[1]
  data <- sim$data[!(sim$data$protein_Id == protein & sim$data$group_ == "Ctrl"), ]
  lfq <- prolfqua::LFQData$new(data, sim$config)
  observed <- data[data$protein_Id == protein, ]
  means <- tapply(observed$abundance, observed$group_, mean)

  for (method in c("lm", "limma", "limma_voom", "deqms_voom")) {
    res <- protein_contrasts(lfq, "~ group_", one_factor_contrasts, method, protein)
    expect_equal(res$diff, c(B_vs_A = unname(means["B"] - means["A"])), info = method)
    expect_equal(res$missing, c("A_vs_Ctrl", "B_vs_Ctrl"), info = method)
  }
})

test_that("a 2x2 protein without its reference cell gets the estimable cell contrasts in every method", {
  sim <- prolfqua::sim_lfq_data_2factor_config(n_proteins = 30, with_missing = FALSE)
  protein <- sort(unique(sim$data$protein_Id))[1]
  data <- sim$data[!(sim$data$protein_Id == protein & sim$data$Treatment == "A" & sim$data$Background == "X"), ]
  lfq <- prolfqua::LFQData$new(data, sim$config)
  observed <- data[data$protein_Id == protein, ]
  cells <- tapply(observed$abundance, paste0(observed$Treatment, observed$Background), mean)
  contrasts <- c(
    background_in_A = "`TreatmentA:BackgroundZ` - `TreatmentA:BackgroundX`",
    background_in_B = "`TreatmentB:BackgroundZ` - `TreatmentB:BackgroundX`",
    treatment_in_Z = "`TreatmentB:BackgroundZ` - `TreatmentA:BackgroundZ`",
    treatment = "TreatmentB - TreatmentA"
  )
  expected <- c(
    background_in_B = unname(cells["BZ"] - cells["BX"]),
    treatment_in_Z = unname(cells["BZ"] - cells["AZ"])
  )

  for (method in c("lm", "limma", "limma_voom", "deqms_voom")) {
    res <- protein_contrasts(lfq, "~ Treatment * Background", contrasts, method, protein)
    expect_setequal(names(res$diff), names(expected))
    expect_equal(res$diff[names(expected)], expected, info = method)
    # AX is empty, and the main effect averages over it
    expect_equal(res$missing, c("background_in_A", "treatment"), info = method)
  }
})

test_that("limma contrasts and their unmoderated standard errors equal lm's for proteins with missing values", {
  sim <- prolfqua::sim_lfq_data_protein_config(n_proteins = 100, weight_missing = 0.3, seed = 2)
  lfq <- prolfqua::LFQData$new(sim$data, sim$config)
  fit <- function(method) {
    fa <- suppressMessages(prolfqua::build_contrast_analysis(lfq, "~ group_", one_factor_contrasts, method = method))
    suppressMessages(suppressWarnings(fa$get_contrasts()))
  }
  res_limma <- fit("limma")
  res_lm <- fit("lm")
  joined <- dplyr::inner_join(res_limma, res_lm, by = c("protein_Id", "contrast"), suffix = c(".limma", ".lm"))
  # both report the same protein and contrast pairs
  expect_equal(nrow(joined), nrow(res_limma))
  expect_equal(nrow(joined), nrow(res_lm))
  # the comparison covers proteins with missing values, where contrasts.fit() standard errors are approximate
  long <- lfq$data_long()
  with_missing <- unique(long$protein_Id[is.na(long[[lfq$response()]])])
  expect_gt(sum(joined$protein_Id %in% with_missing), 20)
  expect_equal(joined$diff.limma, joined$diff.lm)
  expect_equal(unname(joined$std.error.unmoderated.limma), joined$std.error.unmoderated.lm)
  expect_equal(joined$df.unmoderated.limma, joined$df.unmoderated.lm)
})
