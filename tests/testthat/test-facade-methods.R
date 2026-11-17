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

test_that("printing a facade is informative and does not evaluate contrasts", {
  facade <- facade_methods_fit()
  facade$contrast <- list(get_contrasts = function(...) stop("contrasts must remain lazy"))
  expect_output(returned <- withVisible(print(facade)), "prolfqua analysis: lm")
  expect_output(print(facade), "Formula:.*group_")
  expect_output(print(facade), "Input subjects: 12")
  expect_identical(returned$value, facade)
  expect_false(returned$visible)
})

test_that("summary counts agree with the existing contrast interface", {
  facade <- facade_methods_fit()
  result <- summary(facade, FDR_threshold = 0.2, diff_threshold = 0.3)
  expect_s3_class(result, "summary.ContrastsFacadeBase")
  counts <- result$per_contrast
  selected <- facade$filter_significant(FDR_threshold = 0.2, diff_threshold = 0.3)
  expect_equal(counts$estimated, nrow(facade$get_contrasts()))
  expect_equal(counts$missing, nrow(facade$get_missing()))
  expect_equal(counts$significant, nrow(selected))
  expect_equal(counts$up, sum(selected$diff > 0))
  expect_equal(counts$down, sum(selected$diff < 0))
  expect_output(print(result), "prolfqua analysis summary: lm")
})

test_that("hurdle printing and summaries describe both component models", {
  sim <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.4, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  expect_warning(
    facade <- facade_methods_fit("hurdle_nested", lfq),
    "Zero sample variances detected"
  )
  expect_output(print(facade), "intensity formula:")
  expect_output(print(facade), "detection formula:")
  result <- summary(facade)
  expect_equal(result$fit_status$component, c("intensity", "detection"))
  expect_equal(result$per_contrast$estimated, nrow(facade$get_contrasts()))
})
