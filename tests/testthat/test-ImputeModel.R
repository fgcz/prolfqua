# Test imputed_model S3 class and build_model_impute

test_that("coefficient tables use the borrowed covariance and df of imputed fits", {
  data <- transform(iris, protein_Id = "p1")
  model <- build_model(data, strategy_lm("Sepal.Length ~ Species"), subject_id = "protein_Id")
  fit <- model$model_df$linear_model[[1]]
  covariance <- stats::vcov(fit) * 10
  model$model_df$linear_model[[1]] <- new_imputed_model(
    fit,
    covariance,
    borrowed_sigma = 1,
    borrowed_df = 7,
    n_observed = nrow(data)
  )
  actual <- model$get_coefficients()
  expect_equal(actual$Std..Error, unname(sqrt(diag(covariance))))
  expect_equal(actual$df, rep(7, nrow(actual)))
  expect_equal(
    actual$Pr...t..,
    unname(2 * stats::pt(-abs(actual$Estimate / sqrt(diag(covariance))), df = 7))
  )
})

test_that("imputed_model S3 dispatch works", {
  withr::local_seed(1)
  dat <- data.frame(x = 1:10, y = rnorm(10))
  m <- lm(y ~ x, data = dat)

  borrowed_vcov <- matrix(c(1, 0.1, 0.1, 0.5), nrow = 2, dimnames = list(names(coef(m)), names(coef(m))))
  borrowed_sigma <- 0.42
  borrowed_df <- 7

  wrapped <- new_imputed_model(m, borrowed_vcov, borrowed_sigma, borrowed_df, n_observed = 8)

  expect_s3_class(wrapped, "imputed_model")
  expect_s3_class(wrapped, "lm")
  expect_equal(vcov(wrapped), borrowed_vcov)
  expect_equal(sigma(wrapped), borrowed_sigma)
  expect_equal(df.residual(wrapped), borrowed_df)
  expect_equal(coefficients(wrapped), coefficients(m))
})


test_that("build_model_impute recovers singular proteins", {
  istar <- sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0.5, seed = 42)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  response <- lfqdata$response()
  strat <- strategy_lm(paste(response, "~ group_"))

  mod_no_impute <- build_model(lfqdata, strat)
  mod_impute <- build_model_impute(lfqdata, strat)

  n_good_no_impute <- sum(
    mod_no_impute$model_df$has_model_fit &
      !is.na(mod_no_impute$model_df$isSingular) &
      !mod_no_impute$model_df$isSingular,
    na.rm = TRUE
  )
  n_good_impute <- sum(
    mod_impute$model_df$has_model_fit & !is.na(mod_impute$model_df$isSingular) & !mod_impute$model_df$isSingular,
    na.rm = TRUE
  )
  expect_gte(n_good_impute, n_good_no_impute)
  expect_true(grepl("Imputed", mod_impute$model_name))
})


test_that("Contrasts work on imputed model", {
  istar <- sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0.5, seed = 42)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  strat <- strategy_lm(paste(lfqdata$response(), "~ group_"))
  mod <- build_model_impute(lfqdata, strat)

  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  contrastX <- Contrasts$new(mod, Contr)
  res <- contrastX$get_contrasts()

  expect_s3_class(res, "data.frame")
  expect_true(all(c("contrast", "diff", "statistic", "p.value", "FDR") %in% colnames(res)))
  expect_true(all(!is.na(res$statistic)))
})


test_that("LOD-imputed proteins are flagged via estimate_type", {
  muffle_warnings_matching(one_group_warning, expr = {
    # Fixture chosen to guarantee some singular fits that get rescued via
    # impute_refit_singular. Without the rescue, build_model alone would
    # leave those proteins with NA estimates; with build_model_impute,
    # they fit successfully but should be tagged so consumers can tell
    # rescued rows from clean ones.
    istar <- sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.5, seed = 42)
    lfqdata <- LFQData$new(istar$data, istar$config)
    lfqdata$rename_response("transformedIntensity")
    strat <- strategy_lm(paste(lfqdata$response(), "~ group_"))

    mod_plain <- build_model(lfqdata, strat)
    mod_impute <- build_model_impute(lfqdata, strat)

    # build_model_impute should add an imputed column flagging refit rows.
    expect_true("imputed" %in% colnames(mod_impute$model_df))
    n_imputed <- sum(mod_impute$model_df$imputed, na.rm = TRUE)
    expect_gt(n_imputed, 0)

    # And the underlying plain model_df should not have it (back-compat).
    expect_false("imputed" %in% colnames(mod_plain$model_df))

    Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

    # modelName is the model identity (uniform); the LOD-rescue state lives in
    # the separate estimate_type column ("lod_imputed" for refit proteins).
    ctr_impute <- Contrasts$new(mod_impute, Contr)$get_contrasts()
    expect_setequal(unique(ctr_impute$modelName), "WaldTest")
    expect_true("estimate_type" %in% colnames(ctr_impute))
    expect_setequal(unique(ctr_impute$estimate_type), c("observed", "lod_imputed"))
    expect_equal(sum(ctr_impute$estimate_type == "lod_imputed"), n_imputed)

    # The plain Contrasts path: a single uniform modelName, all observed.
    ctr_plain <- Contrasts$new(mod_plain, Contr)$get_contrasts()
    expect_setequal(unique(ctr_plain$modelName), "WaldTest")
    expect_setequal(unique(ctr_plain$estimate_type), "observed")
  })
})


test_that("lm_impute facade surfaces LOD-rescued rows via estimate_type", {
  muffle_warnings_matching(one_group_warning, expr = {
    istar <- sim_lfq_data_protein_config(n_proteins = 50, weight_missing = 0.5, seed = 42)
    lfqdata <- LFQData$new(istar$data, istar$config)
    lfqdata$rename_response("transformedIntensity")
    contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

    fa <- ContrastsLMImputeFacade$new(lfqdata, "~ group_", contrasts)
    res <- fa$get_contrasts()
    # modelName is the facade key; rescued rows are distinguished by estimate_type.
    expect_setequal(unique(res$modelName), "lm_impute")
    expect_true("estimate_type" %in% colnames(res))
    expect_true("observed" %in% res$estimate_type)
    expect_true("lod_imputed" %in% res$estimate_type)
  })
})


test_that("ContrastsLMImputeFacade works end-to-end", {
  istar <- sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0.5, seed = 42)
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
  istar <- sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0.5, seed = 42)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  # some fits of this fixture lack a coefficient, so vcov borrowing falls back to sigma
  expect_warning(
    fa <- ContrastsLMImputeFacade$new(
      lfqdata,
      "~ group_",
      contrasts,
      borrow_method = "vcov"
    ),
    "using the sigma method instead"
  )
  res <- fa$get_contrasts()
  expect_s3_class(res, "data.frame")
  # imputation gives every protein an estimate
  expect_equal(nrow(res), dplyr::n_distinct(lfqdata$data_long()$protein_Id))
})


test_that("compute_borrowed_variance keeps fits missing a whole group as donors", {
  bb <- sim_lfq_data_protein_config(n_proteins = 10, with_missing = FALSE)
  first <- unique(bb$data$protein_Id)[1]
  bb$data$abundance[bb$data$protein_Id == first & bb$data$group_ == "A"] <- NA
  mod <- build_model(bb$data, strategy_lm("abundance ~ group_"), subject_id = "protein_Id")
  md <- mod$model_df
  expect_false(md$isSingular[md$protein_Id == first])
  borrowed <- compute_borrowed_variance(md, method = "sigma")
  expect_equal(borrowed$sigma, stats::median(md$sigma))
})


test_that("df_method borrowed works", {
  istar <- sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0.5, seed = 42)
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
  expect_equal(nrow(res), dplyr::n_distinct(lfqdata$data_long()$protein_Id))
})

test_that("imputed fits identify observed and completed rows", {
  sim <- sim_lfq_data_protein_config(n_proteins = 12, weight_missing = 0, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)
  dat <- lfq$data_long()
  protein <- unique(dat$protein_Id)[1]
  lfq$set_data(dat[!(dat$protein_Id == protein & dat$group_ == "A"), ])
  model <- build_model_impute(lfq, strategy_lm(paste(lfq$response(), "~ group_")))
  row <- which(model$model_df$protein_Id == protein)
  fit <- model$model_df$linear_model[[row]]
  data <- model$model_df$data[[row]]
  expect_equal(attr(fit, "imputed_rows"), data$group_ == "A")
  expect_equal(sum(!attr(fit, "imputed_rows")), attr(fit, "n_observed"))
})

test_that("rfit reports its fitted observation count", {
  skip_if_not_installed("Rfit")
  dat <- iris
  dat$Sepal.Length[1] <- NA_real_
  fit <- strategy_rfit("Sepal.Length ~ Species")$model_fun(dat)
  expect_equal(stats::nobs(fit), 149L)
})

test_that("build_model_impute rejects input and strategies outside its contract", {
  sim <- sim_lfq_data_peptide_config(n_proteins = 4, with_missing = FALSE, seed = 1)
  cfg <- sim$config$clone(deep = TRUE)
  cfg$hierarchy_depth <- 1
  peptide_level <- LFQData$new(sim$data, cfg)
  expect_error(
    build_model_impute(peptide_level, strategy_lm(paste(peptide_level$response(), "~ group_"))),
    "requires aggregated LFQData",
    class = "prolfqua_error"
  )

  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  protein_sim <- sim_lfq_data_protein_config(n_proteins = 4, with_missing = FALSE, seed = 1, paired = TRUE)
  protein_level <- LFQData$new(protein_sim$data, protein_sim$config)
  expect_error(
    build_model_impute(protein_level, strategy_lmer(paste(protein_level$response(), "~ group_ + (1 | subject_)"))),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("build_model_impute returns the imputed column also when no fit needs a refit", {
  istar <- sim_lfq_data_protein_config(n_proteins = 10, weight_missing = 0, seed = 1)
  lfqdata <- LFQData$new(istar$data, istar$config)
  strategy <- strategy_lm(paste(lfqdata$response(), "~ group_"))
  model <- suppressMessages(build_model_impute(lfqdata, strategy))
  expect_true("imputed" %in% colnames(model$model_df))
  expect_false(any(model$model_df$imputed))
})
