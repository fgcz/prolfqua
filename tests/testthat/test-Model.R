# Test Model and ModelFirth R6 classes

test_that("Model (lm strategy)", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20)
  config <- istar$config

  modelFunction <- strategy_lm("abundance ~ group_")
  mod <- build_model(
    istar$data,
    modelFunction,
    subject_id = config$hierarchy_keys_depth()
  )

  # get_coefficients
  coefs <- mod$get_coefficients()
  expect_s3_class(coefs, "data.frame")
  expect_true("Estimate" %in% colnames(coefs))
  expect_true("factor" %in% colnames(coefs))

  # get_anova
  anova <- mod$get_anova()
  expect_s3_class(anova, "data.frame")
  expect_true("p.value" %in% colnames(anova))

  # coef_histogram
  ch <- mod$coef_histogram()
  expect_type(ch, "list")
  expect_s3_class(ch$plot, "ggplot")

  # coef_volcano
  cv <- mod$coef_volcano()
  expect_type(cv, "list")
  expect_s3_class(cv$plot, "ggplot")

  # coef_pairs
  cp <- mod$coef_pairs()
  expect_type(cp, "list")
  expect_s3_class(cp$plot, "data.frame")

  # anova_histogram
  ah <- mod$anova_histogram()
  expect_type(ah, "list")
  expect_s3_class(ah$plot, "ggplot")
})

test_that("strategy_rlm sigma matches the scale embedded in vcov() (moderation coherence)", {
  set.seed(1)
  d <- data.frame(
    abundance = c(stats::rnorm(10), stats::rnorm(10, mean = 2)),
    group_ = rep(c("A", "B"), each = 10)
  )
  strat <- strategy_rlm("abundance ~ group_")
  fit <- strat$model_fun(d)

  # ContrastsModerated rescales the raw Wald statistic by sigma / sqrt(var.post).
  # That is only coherent when the strategy's sigma equals the scale factor
  # embedded in vcov() (i.e. std.error == sigma * unscaled). MASS::rlm builds its
  # vcov()/summary() standard errors from the robust scale fit$s, so the strategy
  # sigma must be fit$s -- not stats::sigma(), which is the larger
  # ordinary-residual scale and would leave a spurious per-protein factor in
  # every moderated rlm p-value.
  se_vcov <- sqrt(diag(vcov(fit)))
  unscaled <- se_vcov / fit$s # scale-free contrast geometry
  expect_equal(strat$sigma(fit) * unscaled, se_vcov)
  expect_equal(strat$sigma(fit), fit$s)
})

test_that("Model (lm strategy with weights)", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20)
  config <- istar$config

  # Add a weight column with non-uniform weights
  set.seed(42)
  istar$data$wt_col <- runif(nrow(istar$data), min = 0.1, max = 10)

  # Fit unweighted model
  mod_unw <- build_model(
    istar$data,
    strategy_lm("abundance ~ group_"),
    subject_id = config$hierarchy_keys_depth()
  )

  # Fit weighted model
  mod_w <- build_model(
    istar$data,
    strategy_lm("abundance ~ group_", weights = "wt_col"),
    subject_id = config$hierarchy_keys_depth()
  )

  coefs_unw <- mod_unw$get_coefficients()
  coefs_w <- mod_w$get_coefficients()

  expect_s3_class(coefs_w, "data.frame")
  expect_true("Estimate" %in% colnames(coefs_w))

  # Weighted and unweighted estimates should differ
  expect_false(
    all(coefs_unw$Estimate == coefs_w$Estimate),
    info = "Weighted and unweighted models should produce different estimates"
  )
})

test_that("ModelFirth", {
  istar <- sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
  lfq <- LFQData$new(istar$data, istar$config)
  lfq$set_data(encode_bin_resp(lfq))
  lfq$set_config_value("bin_resp", "bin_resp")
  formula <- paste0(lfq$get_config()$bin_resp, "~ group_")
  mod <- build_model_logistf(lfq, formula)

  coefs <- mod$get_coefficients()
  expect_s3_class(coefs, "data.frame")
  expect_true("Estimate" %in% colnames(coefs))

  ch <- mod$coef_histogram()
  expect_type(ch, "list")
  expect_s3_class(ch$plot, "ggplot")

  cp <- mod$coef_pairs()
  expect_type(cp, "list")

  cv <- mod$coef_volcano()
  expect_type(cv, "list")
  expect_s3_class(cv$plot, "ggplot")
})

test_that("StrategyLogistf avoids coefficient-wise profile likelihood fits", {
  data <- data.frame(
    bin_resp = c(1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L),
    condition = rep(c("A", "B"), each = 4)
  )

  fit <- strategy_logistf("bin_resp ~ condition")$model_fun(data)

  expect_setequal(fit$method.ci, "Wald")
})

test_that("build_model_glm_protein and build_model_glm_peptide return ModelFirth", {
  prot <- sim_lfq_data_protein_config(Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 7)
  prot_lfq <- LFQData$new(prot$data, prot$config)
  prot_mod <- build_model_glm_protein(prot_lfq, "~ group_")
  expect_true(inherits(prot_mod, "ModelFirth"))
  expect_s3_class(prot_mod$get_coefficients(), "data.frame")

  pep <- sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 7)
  pep_lfq <- LFQData$new(pep$data, pep$config)
  pep_mod <- build_model_glm_peptide(pep_lfq, "~ group_")
  expect_true(inherits(pep_mod, "ModelFirth"))
  expect_s3_class(pep_mod$get_coefficients(), "data.frame")
})
