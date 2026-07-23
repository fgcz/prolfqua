test_that("binomial counts reuse the completed Firth missingness representation", {
  istar <- sim_lfq_data_peptide_config(Nprot = 12, weight_missing = 1, seed = 5)
  lfqdata <- LFQData$new(istar$data, istar$config)
  response <- lfqdata$response()

  sparse <- lfqdata$get_copy()
  sparse$set_data(sparse$data_long()[!is.na(sparse$data_long()[[response]]), ])

  prepared <- .prepare_detection_lfqdata(sparse)
  expect_setequal(unique(prepared$data_long()$bin_resp), c(0L, 1L))

  counts <- .summarize_detection_counts(prepared)
  grouping <- unique(c(
    prepared$subject_id(),
    prepared$sample_name(),
    prepared$isotope_label(),
    prepared$factor_keys()
  ))
  expected <- prepared$data_long() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) |>
    dplyr::summarise(
      detected = sum(.data$bin_resp),
      available = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(undetected = .data$available - .data$detected) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(grouping)))

  expect_equal(
    counts |> dplyr::arrange(dplyr::across(dplyr::all_of(grouping))),
    expected
  )
  expect_true(any(counts$detected == 0L))
  expect_true(all(counts$detected + counts$undetected == counts$available))

  denominator_counts <- counts |>
    dplyr::group_by(dplyr::across(dplyr::all_of(prepared$subject_id()))) |>
    dplyr::summarise(n_available = dplyr::n_distinct(.data$available), .groups = "drop")
  expect_true(all(denominator_counts$n_available == 1L))
})

test_that("strategy_binomial matches an independent quasibinomial fit", {
  dat <- data.frame(
    group_ = factor(rep(c("A", "B"), each = 4)),
    detected = c(1, 2, 1, 3, 4, 5, 3, 5),
    undetected = c(4, 3, 4, 2, 1, 0, 2, 0)
  )

  strategy <- strategy_binomial("~ group_", prior_count = 0.1)
  fit <- strategy$model_fun(dat)
  reference <- stats::glm(
    cbind(detected + 0.1, undetected + 0.1) ~ group_,
    family = stats::quasibinomial(),
    data = dat
  )

  expect_s3_class(fit, "glm")
  expect_equal(stats::coef(fit), stats::coef(reference))
  expect_equal(stats::vcov(fit), stats::vcov(reference))
  expect_equal(strategy$sigma(fit), sqrt(summary(fit)$dispersion))
  expect_equal(strategy$sigma(fit), stats::sigma(fit))
  expect_equal(strategy$df_residual(fit), stats::df.residual(fit))

  unscaled <- sqrt(diag(summary(fit)$cov.unscaled))
  expect_equal(strategy$sigma(fit) * unscaled, sqrt(diag(stats::vcov(fit))))

  linfct <- matrix(
    c(0, 1),
    nrow = 1,
    dimnames = list("B_vs_A", names(stats::coef(fit)))
  )
  contrast <- compute_contrast(fit, linfct)
  expect_equal(contrast$estimate, unname(stats::coef(fit)[["group_B"]]))
  expect_equal(
    contrast$std.error,
    as.numeric(sqrt(linfct %*% stats::vcov(fit) %*% t(linfct)))
  )
  expect_false(strategy$isSingular(fit))
  expect_false(strategy$is_mixed)
  expect_equal(strategy$model_name, "binomial_nested")
  expect_error(strategy_binomial("~ group_", prior_count = Inf), "non-negative number")
})

test_that("strategy_binomial pseudocount bounds separation and supports one child", {
  separated <- data.frame(
    group_ = factor(rep(c("A", "B"), each = 4)),
    detected = c(rep(0, 4), rep(5, 4)),
    undetected = c(rep(5, 4), rep(0, 4))
  )
  fit_small <- strategy_binomial("~ group_", prior_count = 0.1)$model_fun(separated)
  fit_large <- strategy_binomial("~ group_", prior_count = 0.5)$model_fun(separated)

  expect_true(all(is.finite(stats::coef(fit_small))))
  expect_equal(unname(stats::coef(fit_small)[["group_B"]]), 2 * log(51), tolerance = 1e-6)
  expect_gt(abs(stats::coef(fit_small)[["group_B"]]), abs(stats::coef(fit_large)[["group_B"]]))

  one_child <- data.frame(
    group_ = factor(rep(c("A", "B"), each = 4)),
    detected = c(0, 1, 0, 1, 1, 1, 0, 1),
    undetected = c(1, 0, 1, 0, 0, 0, 1, 0)
  )
  fit_one <- strategy_binomial("~ group_", prior_count = 0.1)$model_fun(one_child)
  expect_s3_class(fit_one, "glm")
  expect_equal(stats::df.residual(fit_one), 6)
  expect_true(is.finite(stats::sigma(fit_one)))
})

test_that("strategy_binomial validates its pseudocount", {
  invalid_prior_counts <- list(
    "not numeric",
    c(0.1, 0.2),
    NA_real_,
    Inf,
    -0.1
  )

  for (prior_count in invalid_prior_counts) {
    expect_error(
      strategy_binomial("~ group_", prior_count = prior_count),
      "one non-negative number",
      fixed = TRUE
    )
  }
})

test_that("strategy_binomial rejects rank-deficient and insufficient-df fits", {
  rank_deficient <- data.frame(
    group_ = factor(rep(c("A", "B"), each = 3)),
    batch = factor(rep(c("x", "y"), each = 3)),
    detected = c(1, 2, 1, 4, 3, 4),
    undetected = 5 - c(1, 2, 1, 4, 3, 4)
  )
  expect_warning(
    fit_rank <- strategy_binomial("~ group_ + batch")$model_fun(rank_deficient),
    "rank-deficient"
  )
  expect_type(fit_rank, "character")

  too_small <- data.frame(
    group_ = factor(c("A", "A", "B")),
    detected = c(1, 2, 4),
    undetected = c(4, 3, 1)
  )
  expect_warning(
    fit_small <- strategy_binomial("~ group_")$model_fun(too_small),
    "residual degrees of freedom"
  )
  expect_type(fit_small, "character")

  absent_level <- data.frame(
    group_ = factor(rep("A", 4), levels = c("A", "B")),
    detected = c(1, 2, 3, 2),
    undetected = c(4, 3, 2, 3)
  )
  expect_warning(
    fit_absent <- strategy_binomial("~ group_")$model_fun(absent_level),
    "factors with 2 or more levels"
  )
  expect_type(fit_absent, "character")
})

test_that("moderated_p_limma applies an optional posterior variance floor", {
  contrast_df <- data.frame(
    sigma = c(0.20, 0.25, 0.18, 0.22, 0.24),
    df = rep(6, 5),
    statistic = c(2.1, -1.8, 0.5, 3.0, -2.2),
    diff = c(0.8, -0.6, 0.2, 1.2, -0.9)
  )

  unbounded <- moderated_p_limma(contrast_df)
  bounded <- moderated_p_limma(contrast_df, variance_floor = 1)

  expect_true(any(unbounded$moderated.var.post < 1))
  expect_true(all(bounded$moderated.var.post >= 1))
  bounded_rows <- unbounded$moderated.var.post < 1
  expect_true(all(is.infinite(bounded$moderated.df.total[bounded_rows])))
  expect_equal(
    bounded$moderated.statistic[bounded_rows],
    contrast_df$statistic[bounded_rows] * contrast_df$sigma[bounded_rows]
  )

  unchanged <- moderated_p_limma(contrast_df, variance_floor = NULL)
  expect_equal(unchanged, unbounded)
  expect_error(moderated_p_limma(contrast_df, variance_floor = Inf), "positive number")
})

test_that("binomial_nested facade reuses generic model and contrast interfaces", {
  istar <- sim_lfq_data_peptide_config(Nprot = 30, weight_missing = 0.7, seed = 9)
  lfqdata <- LFQData$new(istar$data, istar$config)
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  facade <- build_contrast_analysis(
    lfqdata,
    "~ group_",
    contrasts,
    method = "binomial_nested"
  )

  expect_r6_class(facade, "ContrastsBinomialNestedFacade")
  expect_r6_class(facade$model, "Model")
  expect_r6_class(facade$model$model_strategy, "StrategyBinomial")
  expect_r6_class(facade$contrast, "ContrastsModerated")
  expect_equal(facade$model$subject_id, lfqdata$subject_id())
  config <- facade$get_config()
  expect_equal(config$effect_col, "diff")
  expect_equal(config$avg_abundance_col, "avgAbd")
  expect_false(config$supports_dea_qc)

  result <- facade$get_contrasts()
  expected_columns <- ContrastsInterface$new()$column_description()$column_name
  expect_true(all(expected_columns %in% colnames(result)))
  expect_setequal(unique(result$modelName), "binomial_nested")
  expect_setequal(unique(result$estimate_type), "observed")
  expect_setequal(unique(result$contrast), names(contrasts))
  expect_false(anyNA(result$diff))

  inner <- facade$contrast$Contrast$get_contrasts()
  squeezed <- limma::squeezeVar(inner$sigma^2, df = inner$df)
  expected_var <- pmax(squeezed$var.post, 1)
  expected_statistic <- inner$statistic * inner$sigma / sqrt(expected_var)
  expect_equal(result$statistic, expected_statistic)
  expect_true(all(result$sigma >= 1))

  expect_s3_class(facade$model$get_coefficients(), "data.frame")
  expect_s3_class(facade$model$get_anova(), "data.frame")
  expect_type(facade$model$coef_histogram(), "list")
  expect_type(facade$model$coef_volcano(), "list")
  expect_type(facade$model$coef_pairs(), "list")
  expect_type(facade$model$anova_histogram(), "list")
  expect_r6_class(facade$get_Plotter(), "ContrastsPlotter")
  expect_s3_class(facade$to_wide(), "data.frame")
  expect_s3_class(facade$get_missing(), "data.frame")

  registry_entry <- lookup_facade("binomial_nested")
  expect_equal(registry_entry$class, "ContrastsBinomialNestedFacade")
  expect_equal(registry_entry$needs, "nested")
})

test_that("binomial_nested supports multifactor contrasts and rejects aggregated input", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.8, seed = 11)
  istar$data <- istar$data |>
    dplyr::mutate(
      Treatment = dplyr::if_else(.data$group_ == "B", "B", "A"),
      Background = dplyr::if_else(grepl("_V[34]$", .data$sampleName), "Z", "X")
    )
  istar$config$factors <- list(Treatment = "Treatment", Background = "Background")
  istar$config$factor_depth <- 2
  lfqdata <- LFQData$new(istar$data, istar$config)
  contrasts <- c(
    "treatment" = "TreatmentB - TreatmentA",
    "background" = "BackgroundZ - BackgroundX",
    "background_in_A" = "`TreatmentA:BackgroundZ` - `TreatmentA:BackgroundX`",
    "background_in_B" = "`TreatmentB:BackgroundZ` - `TreatmentB:BackgroundX`",
    "interaction_effect" = "background_in_B - background_in_A"
  )

  facade <- build_contrast_analysis(
    lfqdata,
    "~ Treatment * Background",
    contrasts,
    method = "binomial_nested",
    binomial_bound = FALSE
  )
  result <- facade$get_contrasts()
  expect_true(nrow(result) > 0)
  expect_setequal(unique(result$contrast), names(contrasts))

  protein <- sim_lfq_data_protein_config(Nprot = 10, seed = 4)
  protein_lfq <- LFQData$new(protein$data, protein$config)
  expect_error(
    ContrastsBinomialNestedFacade$new(protein_lfq, "~ group_", c("A_vs_Ctrl" = "group_A - group_Ctrl")),
    "requires LFQData with additional hierarchy"
  )
})
