# Test Model and ModelFirth R6 classes

test_that("Model (lm strategy)", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20)
  config <- istar$config

  modelFunction <- strategy_lm("abundance ~ group_")
  mod <- build_model(
    istar$data,
    modelFunction,
    subject_Id = config$hierarchy_keys_depth()
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

test_that("ModelFirth", {
  istar <- sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE,
    weight_missing = 0.5, seed = 3)
  istar$data <- encode_bin_resp(istar$data, istar$config)
  lfq <- LFQData$new(istar$data, istar$config)
  formula <- paste0(lfq$config$bin_resp, "~ group_")
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
