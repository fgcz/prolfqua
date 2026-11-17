facade_methods_input <- function() {
  sim <- sim_lfq_data_protein_config(Nprot = 12, weight_missing = 0, seed = 42)
  LFQData$new(sim$data, sim$config)
}

facade_methods_fit <- function(method = "lm", lfq = facade_methods_input()) {
  build_contrast_analysis(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"), method)
}

test_that("internal fit accessors align omitted observations and aliased coefficients", {
  data <- iris
  data$Sepal.Length[c(2, 8)] <- NA_real_
  data$duplicate <- data$Sepal.Width
  fit <- stats::lm(Sepal.Length ~ Sepal.Width + duplicate, data)
  coefficients <- .fit_coef(fit)
  expect_equal(coefficients$estimate, unname(stats::coef(fit)))
  expect_equal(coefficients$std.error, unname(sqrt(diag(stats::vcov(fit)))))
  values <- .fit_values(fit, data)
  expect_equal(values$row, setdiff(seq_len(nrow(data)), c(2L, 8L)))
  expect_equal(values$fitted, unname(stats::fitted(fit)))
  expect_equal(values$response, values$fitted + values$residual)
  expect_equal(.fit_nobs(fit), 148L)
})

test_that("model accessors return ungrouped keyed tables and fit status", {
  facade <- facade_methods_fit()
  model <- facade$model
  coefficients <- .model_coef(model)
  expect_named(coefficients, c("protein_Id", "term", "estimate", "std.error", "df"))
  expect_false(dplyr::is_grouped_df(coefficients))
  expect_equal(nrow(coefficients), 36L)
  expect_equal(.model_nobs(model)$nobs, rep(12, 12))
  values <- .model_values(model, facade$.lfqdata)
  expect_equal(nrow(values), 144L)
  expect_equal(values$fitted + values$residual, values$response)
  expect_false(any(values$imputed))
  expect_equal(.model_fit_status(model)$fitted, 12L)
  expect_equal(.model_formula(model), model$model_strategy$formula)
})

test_that("matrix accessors preserve subject and sample alignment", {
  facade <- facade_methods_fit("limma")
  model <- facade$model
  coefficients <- .model_coef(model)
  expect_equal(coefficients$estimate, as.numeric(model$fit$coefficients))
  expect_equal(.model_nobs(model)$nobs, unname(rowSums(model$observed)))
  values <- .model_values(model, facade$.lfqdata)
  expect_equal(values$fitted, as.numeric(stats::fitted(model$fit)))
  expect_equal(values$response, as.numeric(model$E))
  expect_equal(values$fitted + values$residual, values$response)
})

test_that("imputed accessors use borrowed uncertainty and omit filled responses", {
  lfq <- facade_methods_input()
  data <- lfq$data_long()
  protein <- unique(data$protein_Id)[1]
  lfq$set_data(data[!(data$protein_Id == protein & data$group_ == "A"), ])
  facade <- facade_methods_fit("lm_impute", lfq)
  model <- facade$model
  fit <- model$model_df$linear_model[[which(model$model_df$protein_Id == protein)]]
  coefficients <- .model_coef(model)
  expect_equal(coefficients$std.error[coefficients$protein_Id == protein], unname(sqrt(diag(stats::vcov(fit)))))
  values <- .model_values(model, facade$.lfqdata)
  rescued <- values[values$protein_Id == protein, ]
  expect_true(all(rescued$imputed))
  expect_false(any(rescued$group_ == "A"))
  expect_equal(nrow(rescued), attr(fit, "n_observed"))
})

test_that("failed fits are skipped without losing table schemas", {
  facade <- facade_methods_fit()
  model <- facade$model
  model$model_df$has_model_fit[] <- FALSE
  model$model_df$linear_model <- rep(list("failed"), nrow(model$model_df))
  expect_equal(nrow(.model_coef(model)), 0L)
  expect_equal(nrow(.model_nobs(model)), 0L)
  expect_equal(nrow(.model_values(model, facade$.lfqdata)), 0L)
  expect_equal(.model_fit_status(model)$failed, 12L)
  expect_error(.model_coef(NULL), "x\\$contrast\\$intensity", class = "prolfqua_error_not_available")
})
