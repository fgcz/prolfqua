# Test lm_imputed S3 class and build_model_impute

test_that("lm_imputed S3 dispatch works", {
  dat <- data.frame(x = 1:10, y = rnorm(10))
  m <- lm(y ~ x, data = dat)

  borrowed_vcov <- matrix(c(1, 0.1, 0.1, 0.5), nrow = 2, dimnames = list(names(coef(m)), names(coef(m))))
  borrowed_sigma <- 0.42
  borrowed_df <- 7

  wrapped <- new_lm_imputed(m, borrowed_vcov, borrowed_sigma, borrowed_df, n_observed = 8)

  expect_s3_class(wrapped, "lm_imputed")
  expect_s3_class(wrapped, "lm")
  expect_equal(vcov(wrapped), borrowed_vcov)
  expect_equal(sigma(wrapped), borrowed_sigma)
  expect_equal(df.residual(wrapped), borrowed_df)
  expect_equal(coefficients(wrapped), coefficients(m))
})


test_that("build_model_impute recovers singular proteins", {
  istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5, seed = 42)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  response <- lfqdata$config$get_response()
  strat <- strategy_lm(paste(response, "~ group_"))

  mod_no_impute <- build_model(lfqdata, strat)
  mod_impute <- build_model_impute(lfqdata, strat)

  n_good_no_impute <- sum(
    mod_no_impute$modelDF$exists_lmer & !is.na(mod_no_impute$modelDF$isSingular) & !mod_no_impute$modelDF$isSingular,
    na.rm = TRUE
  )
  n_good_impute <- sum(
    mod_impute$modelDF$exists_lmer & !is.na(mod_impute$modelDF$isSingular) & !mod_impute$modelDF$isSingular,
    na.rm = TRUE
  )
  expect_gte(n_good_impute, n_good_no_impute)
  expect_true(grepl("Imputed", mod_impute$modelName))
})


test_that("Contrasts work on imputed model", {
  istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5, seed = 42)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  strat <- strategy_lm(paste(lfqdata$config$get_response(), "~ group_"))
  mod <- build_model_impute(lfqdata, strat)

  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contrastX <- Contrasts$new(mod, Contr)
  res <- contrastX$get_contrasts()

  expect_s3_class(res, "data.frame")
  expect_true(all(c("contrast", "diff", "statistic", "p.value", "FDR") %in% colnames(res)))
  expect_true(all(!is.na(res$statistic)))
})


test_that("ContrastsLMImputeFacade works end-to-end", {
  istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5, seed = 42)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa_impute <- ContrastsLMImputeFacade$new(lfqdata, "~ group_", contrasts)
  fa_lm <- ContrastsLMFacade$new(lfqdata, "~ group_", contrasts)

  res_impute <- fa_impute$get_contrasts()
  res_lm <- fa_lm$get_contrasts()

  expect_s3_class(res_impute, "data.frame")
  expect_gte(nrow(res_impute), nrow(res_lm))

  missing_impute <- fa_impute$get_missing()
  missing_lm <- fa_lm$get_missing()
  expect_lte(nrow(missing_impute), nrow(missing_lm))

  expect_s3_class(fa_impute$get_Plotter(), "ContrastsPlotter")
  tw <- fa_impute$to_wide()
  expect_s3_class(tw, "data.frame")
})


test_that("borrow_method vcov works", {
  istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5, seed = 42)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- ContrastsLMImputeFacade$new(
    lfqdata,
    "~ group_",
    contrasts,
    borrow_method = "vcov"
  )
  res <- fa$get_contrasts()
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
})


test_that("df_method borrowed works", {
  istar <- sim_lfq_data_protein_config(Nprot = 30, weight_missing = 0.5, seed = 42)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  fa <- ContrastsLMImputeFacade$new(
    lfqdata,
    "~ group_",
    contrasts,
    df_method = "borrowed"
  )
  res <- fa$get_contrasts()
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
})
