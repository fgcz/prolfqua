facade_methods_input <- function() {
  sim <- sim_lfq_data_protein_config(n_proteins = 12, weight_missing = 0, seed = 42)
  LFQData$new(sim$data, sim$config)
}

facade_methods_fit <- function(method = "lm", lfq = facade_methods_input()) {
  build_contrast_analysis(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"), method)
}

# limma::squeezeVar() warns about zero sample variances in limma 3.68 but not in 3.99;
# these tests are about printing and accessors, not about that warning.
muffle_limma_zero_variance <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Zero sample variances detected", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
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
  expect_output(print(facade), "Coefficient preview:")
  expect_output(print(facade), "estimate")
  expect_identical(returned$value, facade)
  expect_false(returned$visible)
})

test_that("printing extracts coefficients from only the first few fits", {
  facade <- facade_methods_fit()
  extract <- .fit_coef
  extracted <- 0L
  local_mocked_bindings(.fit_coef = function(fit) {
    extracted <<- extracted + 1L
    extract(fit)
  })
  expect_output(print(facade), "Coefficient preview:")
  expect_lte(extracted, 6L)
})

test_that("matrix fit status does not infer singularity from failed fits", {
  facade <- facade_methods_fit("limma")
  expect_identical(.model_fit_status(facade$model)$singular, NA_integer_)
  facade$model$fit$coefficients[1, 1] <- NA_real_
  status <- .model_fit_status(facade$model)
  expect_equal(status$fitted, 11L)
  expect_equal(status$failed, 1L)
  expect_identical(status$singular, NA_integer_)
})

test_that("coefficient previews preserve full-table values and ordering", {
  for (method in c("lm", "limma")) {
    facade <- facade_methods_fit(method)
    coefficients <- stats::coef(facade)
    for (n in c(1L, 6L, 15L, 100L)) {
      expect_equal(.model_coef(facade$model, n), utils::head(coefficients, n))
    }
  }
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
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.4, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  facade <- muffle_limma_zero_variance(facade_methods_fit("hurdle_nested", lfq))
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
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.4, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  facade <- muffle_limma_zero_variance(facade_methods_fit("hurdle_nested", lfq))
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
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.4, seed = 42)
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

# These generated data deliberately include single-peptide proteins, sparse
# groups, tied peptide counts and fully detected proteins. Those backend warnings
# are exercised elsewhere; unexpected accessor warnings must still fail here.
facade_fixture_warnings <- paste(
  single_peptide_warning,
  one_group_warning,
  rlm_convergence_warning,
  loess_ties_warning,
  "Zero sample variances detected",
  sep = "|"
)

facade_method_case <- function(method) {
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, seed = 42)
  peptide <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  if (method == "limpa") {
    input <- AggregateLimpa$new(peptide, "protein")$aggregate()
  } else if (lookup_facade(method)$needs == "nested") {
    input <- peptide
  } else {
    input <- peptide$get_Aggregator()$aggregate()
  }
  model <- if (method == "lmer_nested") "~ group_ + (1 | peptide_Id) + (1 | sampleName)" else "~ group_"
  build_contrast_analysis(input, model, c(A_vs_Ctrl = "group_A - group_Ctrl"), method)
}

facade_method_dependencies <- list(
  rfit = "Rfit",
  rfit_impute = "Rfit",
  lmer_nested = c("lme4", "lmerTest"),
  firth = "logistf",
  limpa = "limpa",
  limpa_nested = "limpa"
)

#' @srrstats {RE7.3} The registry loop tests supported accessors, covariance consistency and documented hurdle errors.
#' @noRd
NULL

for (method in names(list_facades())) {
  test_that(paste("standard analysis methods work for", method), {
    for (package in facade_method_dependencies[[method]]) {
      skip_if_not_installed(package)
    }
    muffle_warnings_matching(pattern = facade_fixture_warnings, expr = {
      facade <- facade_method_case(method)
      expect_output(print(facade), paste("prolfqua analysis:", method))
      summary <- summary(facade)
      expect_true(all(c("subjects", "fitted", "failed", "singular", "imputed") %in% names(summary$fit_status)))
      expect_equal(summary$fit_status$fitted + summary$fit_status$failed, summary$fit_status$subjects)
      expect_equal(summary$per_contrast$estimated, nrow(facade$get_contrasts()))
      expect_equal(summary$per_contrast$missing, nrow(facade$get_missing()))
      expect_s3_class(graphics::plot(facade), "ggplot")
      if (method == "hurdle_nested") {
        for (accessor in list(
          stats::coef,
          stats::vcov,
          stats::confint,
          stats::formula,
          stats::nobs,
          stats::fitted,
          stats::residuals
        )) {
          expect_error(accessor(facade), class = "prolfqua_error_not_available")
        }
      } else {
        coefficients <- stats::coef(facade)
        expect_named(coefficients, c(facade$model$subject_id, "term", "estimate", "std.error", "df"))
        expect_gt(nrow(coefficients), 0L)
        covariance <- stats::vcov(facade)
        expect_named(covariance, c(facade$model$subject_id, "term", "term2", "covariance"))
        diagonal <- covariance[covariance$term == covariance$term2, ]
        compared <- dplyr::inner_join(diagonal, coefficients, by = c(facade$model$subject_id, "term"))
        expect_equal(nrow(compared), sum(is.finite(coefficients$estimate)))
        expect_equal(sqrt(compared$covariance), compared$std.error)
        intervals <- stats::confint(facade)
        expect_equal(intervals[, names(coefficients)], coefficients)
        expect_true(all(c("conf.low", "conf.high") %in% names(intervals)))
        expect_s3_class(stats::formula(facade), "formula")
        observations <- stats::nobs(facade)
        expect_named(observations, c(facade$model$subject_id, "nobs"))
        expect_true(all(observations$nobs >= 0))
        values <- stats::fitted(facade)
        residuals <- stats::residuals(facade)
        expect_true(all(
          c(facade$model$subject_id, facade$.lfqdata$sample_name(), "response", "fitted", "imputed") %in% names(values)
        ))
        expect_equal(values$response, values$fitted + residuals$residual)
        expect_false(anyNA(values[[facade$.lfqdata$sample_name()]]))
        if (method %in% c("firth", "binomial_nested", "binomial_firth_nested")) {
          expect_true(all(values$fitted >= 0 & values$fitted <= 1))
        }
        if (method == "lmer_nested") {
          expect_true("peptide_Id" %in% names(values))
          expect_match(paste(deparse(stats::formula(facade)), collapse = ""), "peptide_Id")
          fits <- .model_fits(facade$model)
          for (i in seq_len(nrow(fits))) {
            expected <- lme4::fixef(fits$linear_model[[i]])
            selected <- coefficients[coefficients$protein_Id == fits$protein_Id[i], ]
            expect_equal(selected$term, names(expected))
            expect_equal(selected$estimate, unname(expected))
            table <- stats::coef(summary(fits$linear_model[[i]]))
            expect_equal(selected$std.error, unname(table[, "Std. Error"]))
            expect_equal(selected$df, unname(table[, "df"]))
          }
        }
      }
    })
  })
}

test_that("analysis covariance matches each ordinary backend fit", {
  facade <- facade_methods_fit()
  covariance <- stats::vcov(facade)
  fits <- .model_fits(facade$model)
  for (i in seq_len(nrow(fits))) {
    actual <- covariance[covariance$protein_Id == fits$protein_Id[i], ]
    expected <- stats::vcov(fits$linear_model[[i]])
    expect_equal(actual$covariance, as.numeric(expected))
    expect_equal(actual$term, rep(rownames(expected), ncol(expected)))
    expect_equal(actual$term2, rep(colnames(expected), each = nrow(expected)))
  }
  facade$model$model_df$has_model_fit[] <- FALSE
  expect_named(stats::vcov(facade), names(covariance))
  expect_equal(nrow(stats::vcov(facade)), 0L)
})

test_that("matrix covariance uses each subject's measured cells and precision weights", {
  for (method in c("limma", "limma_voom")) {
    lfq <- facade_methods_input()
    data <- lfq$data_long()
    data[[lfq$response()]][c(1, 14, 27)] <- NA_real_
    lfq$set_data(data)
    facade <- facade_methods_fit(method, lfq)
    model <- facade$model
    covariance <- stats::vcov(facade)
    moderated <- limma::eBayes(model$fit, trend = model$trend, robust = model$robust)
    for (i in seq_len(nrow(model$E))) {
      used <- is.finite(model$E[i, ])
      weights <- if (is.null(model$weights)) rep(1, sum(used)) else model$weights[i, used]
      design <- model$fit$design[used, , drop = FALSE]
      expected <- moderated$s2.post[i] * solve(t(design) %*% diag(weights) %*% design)
      actual <- covariance[covariance$protein_Id == model$rowdata$protein_Id[i], ]
      expect_equal(actual$covariance, as.numeric(expected))
    }
  }
})

test_that("limpa imputation flags describe observations rather than the model label", {
  skip_if_not_installed("limpa")
  facade <- facade_method_case("limpa")
  values <- stats::fitted(facade)
  expected <- rowSums(!facade$model$observed) > 0
  expect_true(any(expected))
  expect_equal(values$imputed, unname(expected[match(values$protein_Id, facade$model$rowdata$protein_Id)]))
  facade$model$model_name <- "custom model label"
  expect_identical(stats::fitted(facade)$imputed, values$imputed)
})

test_that("fitted lm cell means and intervals match prediction for one and two factors", {
  for (two_factors in c(FALSE, TRUE)) {
    lfq <- facade_methods_input()
    formula <- "~ group_"
    if (two_factors) {
      data <- lfq$data_long()
      samples <- data[[lfq$sample_name()]]
      data$batch <- c("X", "Y")[1 + match(samples, unique(samples)) %% 2]
      config <- lfq$get_config()
      config$factors <- c(config$factors, batch = "batch")
      config$factor_depth <- 2
      lfq <- LFQData$new(data, config)
      formula <- "~ group_ * batch"
    }
    facade <- build_contrast_analysis(lfq, formula, c(A_vs_Ctrl = "group_A - group_Ctrl"))
    fits <- .model_fits(facade$model)
    selected <- fits$protein_Id[1:2]
    actual <- .model_group_means(facade$model, selected, level = 0.9)
    for (protein in selected) {
      fit <- fits$linear_model[[match(protein, fits$protein_Id)]]
      factors <- intersect(attr(stats::terms(fit), "term.labels"), names(fit$model))
      data <- unique(fit$model[, factors, drop = FALSE])
      groups <- as.character(.make_interaction_column(data, factors, sep = ":")$interaction)
      expected <- stats::predict(fit, newdata = data, interval = "confidence", level = 0.9)
      rows <- actual[actual$protein_Id == protein, ]
      expected <- expected[match(rows$group, groups), , drop = FALSE]
      expect_equal(rows$estimate, unname(expected[, "fit"]))
      expect_equal(rows$conf.low, unname(expected[, "lwr"]))
      expect_equal(rows$conf.high, unname(expected[, "upr"]))
    }
    plot <- graphics::plot(facade, subject = selected, level = 0.9)
    expect_setequal(plot$data$group, plot$layers[[2]]$data$group)
    expect_setequal(plot$data$protein_Id, selected)
    expect_no_error(ggplot2::ggplot_build(plot))
  }
})

test_that("fitted mixed-model means retain Satterthwaite confidence limits", {
  skip_if_not_installed("lmerTest")
  facade <- muffle_warnings_matching(facade_method_case("lmer_nested"), facade_fixture_warnings)
  fits <- .model_fits(facade$model)
  protein <- fits$protein_Id[1]
  fit <- fits$linear_model[[1]]
  linear <- .linfct_from_model(fit)$linfct_factors
  expected <- .compute_lmer_contrast(fit, linear)
  actual <- .model_group_means(facade$model, protein, level = 0.9)
  margin <- stats::qt(0.95, expected$df) * expected$std.error
  expect_equal(actual$estimate, expected$estimate)
  expect_equal(actual$conf.low, expected$estimate - margin)
  expect_equal(actual$conf.high, expected$estimate + margin)
})

test_that("fitted matrix-model means match moderated contrast intervals", {
  facade <- facade_methods_fit("limma")
  model <- facade$model
  linear <- .linfct_from_model(model$dummy_model)$linfct_factors
  fit <- limma::eBayes(limma::contrasts.fit(model$fit, t(linear)), trend = model$trend, robust = model$robust)
  actual <- .model_group_means(model, model$rowdata$protein_Id, level = 0.95)
  for (i in seq_len(nrow(linear))) {
    expected <- limma::topTable(fit, coef = i, number = Inf, sort.by = "none", confint = TRUE)
    rows <- actual[actual$group == rownames(linear)[i], ]
    expect_equal(rows$estimate, expected$logFC)
    expect_equal(rows$conf.low, expected$CI.L)
    expect_equal(rows$conf.high, expected$CI.R)
  }
})

for (method in c("lm", "limma", "lmer_nested", "firth", "binomial_nested", "binomial_firth_nested", "ropeca_nested")) {
  test_that(paste("fitted plotting renders observations and intervals for", method), {
    for (package in facade_method_dependencies[[method]]) {
      skip_if_not_installed(package)
    }
    facade <- muffle_warnings_matching(facade_method_case(method), facade_fixture_warnings)
    plot <- muffle_warnings_matching(graphics::plot(facade, type = "fitted"), "Zero sample variances detected")
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_true(all(c("response", "subject", "group") %in% names(plot$data)))
    means <- plot$layers[[2]]$data
    expect_true(all(c("estimate", "conf.low", "conf.high", "subject", "group") %in% names(means)))
    if (method %in% c("firth", "binomial_nested", "binomial_firth_nested")) {
      expect_true(all(as.matrix(means[c("estimate", "conf.low", "conf.high")]) >= 0))
      expect_true(all(as.matrix(means[c("estimate", "conf.low", "conf.high")]) <= 1))
    }
    if (method == "ropeca_nested") {
      expect_equal(length(unique(plot$data$subject)), nrow(unique(plot$data[c("protein_Id", "peptide_Id")])))
    }
  })
}

test_that("fitted plots select first-contrast FDRs and validate explicit requests", {
  facade <- facade_methods_fit()
  contrasts <- facade$get_contrasts()
  expected <- head(contrasts$protein_Id[order(contrasts$FDR)], 6)
  expect_setequal(graphics::plot(facade)$data$protein_Id, expected)
  for (subject in list("unknown", NA_character_, character(), 1)) {
    expect_error(graphics::plot(facade, subject = subject), "subject", class = "prolfqua_error_bad_argument")
  }
  for (level in list(0, 1, NA_real_, c(0.9, 0.95))) {
    expect_error(graphics::plot(facade, level = level), "level", class = "prolfqua_error_bad_argument")
  }
  expect_error(graphics::plot(facade, type = c("fitted", "volcano", "ma")), class = "prolfqua_error_bad_argument")
  hurdle <- muffle_warnings_matching(facade_method_case("hurdle_nested"), facade_fixture_warnings)
  expect_equal(graphics::plot(hurdle)$labels, graphics::plot(hurdle, type = "volcano")$labels)
  expect_error(
    graphics::plot(hurdle, type = "fitted"),
    "x\\$contrast\\$intensity",
    class = "prolfqua_error_not_available"
  )
})
