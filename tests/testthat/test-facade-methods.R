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

test_that("coefficient and confidence-interval methods agree with each lm fit", {
  facade <- facade_methods_fit()
  coefficients <- stats::coef(facade)
  intervals <- stats::confint(facade)
  counts <- stats::nobs(facade)
  for (i in seq_len(nrow(facade$model$model_df))) {
    row <- facade$model$model_df[i, ]
    fit <- row$linear_model[[1]]
    selected <- coefficients[coefficients$protein_Id == row$protein_Id, ]
    ci <- intervals[intervals$protein_Id == row$protein_Id, ]
    expect_equal(selected$estimate, unname(stats::coef(fit)))
    expect_equal(unname(as.matrix(ci[, c("conf.low", "conf.high")])), unname(stats::confint(fit)))
    expect_equal(counts$nobs[counts$protein_Id == row$protein_Id], stats::nobs(fit))
  }
  expect_equal(stats::formula(facade), facade$model$model_strategy$formula)
  expect_equal(stats::confint(facade, parm = "group_B"), stats::confint(facade, parm = 2))
  for (level in list(0, 1, NA_real_, Inf, "0.95", c(0.9, 0.95))) {
    expect_error(stats::confint(facade, level = level), class = "prolfqua_error_bad_argument")
  }
  expect_error(stats::confint(facade, parm = "unknown"), class = "prolfqua_error_bad_argument")
  expect_error(stats::confint(facade, parm = -1), class = "prolfqua_error_bad_argument")
})

test_that("public coefficient accessors preserve limma estimates and rfit counts", {
  facade <- facade_methods_fit("limma")
  expect_equal(stats::coef(facade)$estimate, as.numeric(facade$model$fit$coefficients))
  skip_if_not_installed("Rfit")
  facade <- facade_methods_fit("rfit")
  expect_equal(stats::nobs(facade)$nobs, rep(12, 12))
})

test_that("a hurdle result has no single coefficient table or formula", {
  sim <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.4, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  expect_warning(
    facade <- facade_methods_fit("hurdle_nested", lfq),
    "Zero sample variances detected"
  )
  for (method in list(stats::coef, stats::confint, stats::formula, stats::nobs)) {
    expect_error(method(facade), "x\\$contrast\\$detection", class = "prolfqua_error_not_available")
  }
})

test_that("fitted values and response residuals agree with each lm fit", {
  facade <- facade_methods_fit()
  values <- stats::fitted(facade)
  residuals <- stats::residuals(facade)
  expect_equal(values$fitted + residuals$residual, values$response)
  expect_named(values, c("protein_Id", "sampleName", "group_", "response", "fitted", "imputed"))
  for (i in seq_len(nrow(facade$model$model_df))) {
    row <- facade$model$model_df[i, ]
    expect_equal(values$fitted[values$protein_Id == row$protein_Id], unname(stats::fitted(row$linear_model[[1]])))
  }
})

test_that("limma fitted values and residuals use the retained response matrix", {
  facade <- facade_methods_fit("limma")
  expect_equal(stats::fitted(facade)$fitted, as.numeric(stats::fitted(facade$model$fit)))
  expect_equal(stats::residuals(facade)$residual, as.numeric(facade$model$E - stats::fitted(facade$model$fit)))
})

test_that("Firth fitted values are probabilities and residuals are on the response scale", {
  skip_if_not_installed("logistf")
  facade <- facade_methods_fit("firth")
  values <- stats::fitted(facade)
  residuals <- stats::residuals(facade)
  expect_true(all(values$fitted >= 0 & values$fitted <= 1))
  expect_equal(values$response, values$fitted + residuals$residual)
  expect_true(all(values$response %in% c(0, 1)))
})

test_that("public imputed residuals exclude the filled observations", {
  lfq <- facade_methods_input()
  data <- lfq$data_long()
  protein <- unique(data$protein_Id)[1]
  lfq$set_data(data[!(data$protein_Id == protein & data$group_ == "A"), ])
  facade <- facade_methods_fit("lm_impute", lfq)
  values <- stats::fitted(facade)
  expect_false(any(values$protein_Id == protein & values$group_ == "A"))
  expect_true(all(values$imputed[values$protein_Id == protein]))
})

test_that("matrix residuals remain aligned when observed cells are missing", {
  lfq <- facade_methods_input()
  data <- lfq$data_long()
  data[[lfq$response()]][c(1, 14, 27)] <- NA_real_
  lfq$set_data(data)
  facade <- facade_methods_fit("limma", lfq)
  values <- stats::fitted(facade)
  residuals <- stats::residuals(facade)
  cells <- which(facade$model$observed, arr.ind = TRUE)
  expected <- (facade$model$E - stats::fitted(facade$model$fit))[cells]
  expect_equal(residuals$residual, expected)
  expect_equal(values$fitted + residuals$residual, values$response)
  expect_equal(nrow(values), sum(facade$model$observed))
})

test_that("analysis plots return visible ggplots across model families", {
  sim <- sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 0.4, seed = 42)
  peptide <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  for (method in c("lm", "limma", "ropeca_nested", "hurdle_nested")) {
    input <- if (lookup_facade(method)$needs == "nested") peptide else facade_methods_input()
    # The hurdle fixture includes perfectly detected proteins, whose detection
    # sample variances are zero. No other warnings are expected here.
    facade <- muffle_warnings_matching(
      facade_methods_fit(method, input),
      pattern = "Zero sample variances detected"
    )
    for (type in c("volcano", "ma")) {
      output <- withVisible(graphics::plot(facade, type = type))
      expect_s3_class(output$value, "ggplot")
      expect_true(output$visible)
    }
  }
  expect_error(graphics::plot(facade, type = "histogram"), class = "prolfqua_error_bad_argument")
})
