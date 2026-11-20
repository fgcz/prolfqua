#' Coefficient estimates and uncertainty from one successful fit
#' @param fit a fitted backend object
#' @noRd
.fit_coef <- function(fit) UseMethod(".fit_coef")

#' @rdname dot-fit_coef
#' @noRd
#' @export
.fit_coef.default <- function(fit) {
  estimates <- stats::coef(fit)
  errors <- sqrt(diag(stats::vcov(fit)))
  tibble::tibble(
    term = names(estimates),
    estimate = unname(estimates),
    std.error = unname(errors[match(names(estimates), names(errors))]),
    df = stats::df.residual(fit)
  )
}

#' @rdname dot-fit_coef
#' @noRd
#' @export
.fit_coef.merMod <- function(fit) {
  coefficients <- stats::coef(summary(fit))
  estimates <- lme4::fixef(fit)
  tibble::tibble(
    term = names(estimates),
    estimate = unname(estimates),
    std.error = unname(coefficients[names(estimates), "Std. Error"]),
    df = unname(coefficients[names(estimates), "df"])
  )
}

#' @rdname dot-fit_coef
#' @noRd
#' @export
.fit_coef.logistf <- function(fit) {
  estimates <- stats::coef(fit)
  tibble::tibble(
    term = names(estimates),
    estimate = unname(estimates),
    std.error = unname(sqrt(diag(stats::vcov(fit)))[names(estimates)]),
    df = Inf
  )
}

#' Count observations before any LOD refit
#' @param fit a successful fitted backend object
#' @noRd
.fit_nobs <- function(fit) UseMethod(".fit_nobs")

#' @rdname dot-fit_nobs
#' @noRd
#' @export
.fit_nobs.default <- function(fit) stats::nobs(fit)

#' @rdname dot-fit_nobs
#' @noRd
#' @export
.fit_nobs.imputed_model <- function(fit) attr(fit, "n_observed")

#' Convergence of one successful fit
#'
#' @param fit a successful fitted backend object
#' @return list with `converged` (logical; TRUE for non-iterative fits, NA when the backend keeps no status) and
#'   `iterations` (integer; NA for non-iterative fits)
#' @noRd
.fit_convergence <- function(fit) UseMethod(".fit_convergence")

#' @rdname dot-fit_convergence
#' @noRd
#' @export
.fit_convergence.default <- function(fit) {
  if (isS4(fit) && methods::is(fit, "merMod")) {
    conv <- fit@optinfo$conv
    return(list(
      # lme4 sets `code` for convergence problems; a singular fit only adds a message (recorded as singular)
      converged = identical(as.integer(conv$opt), 0L) && all(conv$lme4$code %||% 0L == 0L),
      iterations = as.integer(fit@optinfo$feval)
    ))
  }
  list(converged = TRUE, iterations = NA_integer_)
}

#' @rdname dot-fit_convergence
#' @noRd
#' @export
.fit_convergence.rlm <- function(fit) list(converged = isTRUE(fit$converged), iterations = length(fit$conv))

#' @rdname dot-fit_convergence
#' @noRd
#' @export
.fit_convergence.glm <- function(fit) {
  list(converged = isTRUE(fit$converged), iterations = as.integer(fit$firth_iterations %||% fit$iter))
}

#' @rdname dot-fit_convergence
#' @noRd
#' @export
.fit_convergence.logistf <- function(fit) {
  iterations <- as.integer(fit$iter[[1]])
  list(converged = iterations < fit$control$maxit, iterations = iterations)
}

#' @rdname dot-fit_convergence
#' @noRd
#' @export
.fit_convergence.rfit_prolfqua <- function(fit) list(converged = NA, iterations = NA_integer_)

#' @rdname dot-fit_convergence
#' @noRd
#' @export
.fit_convergence.imputed_model <- function(fit) NextMethod()

#' Fitted responses and response residuals indexed into a subject's data
#' @param fit a successful fitted backend object
#' @param data data frame used for this subject's fit
#' @noRd
.fit_values <- function(fit, data) UseMethod(".fit_values")

#' @rdname dot-fit_values
#' @noRd
#' @export
.fit_values.default <- function(fit, data) {
  fitted <- stats::fitted(fit)
  response <- as.numeric(stats::model.response(stats::model.frame(fit)))
  tibble::tibble(
    row = match(names(fitted), rownames(data)),
    response = response,
    fitted = as.numeric(fitted),
    residual = response - as.numeric(fitted)
  )
}

#' @rdname dot-fit_values
#' @noRd
#' @export
.fit_values.glm <- function(fit, data) {
  fitted <- stats::fitted(fit)
  residual <- stats::residuals(fit, type = "response")
  tibble::tibble(
    row = match(names(fitted), rownames(data)),
    response = as.numeric(fitted + residual),
    fitted = as.numeric(fitted),
    residual = as.numeric(residual)
  )
}

#' @rdname dot-fit_values
#' @noRd
#' @export
.fit_values.rfit_prolfqua <- function(fit, data) {
  fitted <- as.numeric(stats::fitted(fit))
  response <- as.numeric(stats::model.response(fit$model))
  tibble::tibble(
    row = match(rownames(fit$model), rownames(data)),
    response = response,
    fitted = fitted,
    residual = response - fitted
  )
}

#' @rdname dot-fit_values
#' @noRd
#' @export
.fit_values.logistf <- function(fit, data) {
  response <- data[[all.vars(stats::formula(fit))[1]]]
  rows <- which(!is.na(response))
  fitted <- as.numeric(stats::predict(fit, newdata = data[rows, , drop = FALSE], type = "response"))
  tibble::tibble(row = rows, response = response[rows], fitted = fitted, residual = response[rows] - fitted)
}

#' @rdname dot-fit_values
#' @noRd
#' @export
.fit_values.imputed_model <- function(fit, data) {
  values <- NextMethod()
  values[!attr(fit, "imputed_rows")[values$row], , drop = FALSE]
}

#' Successful subject fits without retained dplyr grouping
#' @param model a model adapter
#' @noRd
.model_fits <- function(model) UseMethod(".model_fits")

#' @rdname dot-model_fits
#' @noRd
#' @export
.model_fits.Model <- function(model) {
  dplyr::ungroup(model$model_df[model$model_df$has_model_fit, , drop = FALSE])
}

#' @rdname dot-model_fits
#' @noRd
#' @export
.model_fits.ModelFirth <- function(model) {
  dplyr::ungroup(model$models$model_df[model$models$model_df$has_model_fit, , drop = FALSE])
}

#' Explain why a combined hurdle result has no single fitted model
#' @noRd
.model_not_available <- function() {
  abort_prolfqua(
    "A hurdle analysis has two models; use x$contrast$intensity or x$contrast$detection.",
    class = "prolfqua_error_not_available"
  )
}

#' Coefficient table across model subjects
#' @param model a model adapter
#' @param n maximum number of coefficient rows; finite previews inspect at most n subject fits
#' @noRd
.model_coef <- function(model, n = Inf) UseMethod(".model_coef")

#' Repeat subject identifiers alongside variable-length result tables
#' @param fits ungrouped table with one row per subject
#' @param subject_id names of subject identifier columns
#' @param results list of result tables in subject order
#' @param prototype empty table defining result columns when no fits succeeded
#' @noRd
.bind_fit_results <- function(fits, subject_id, results, prototype) {
  ids <- fits[rep(seq_len(nrow(fits)), vapply(results, nrow, integer(1))), subject_id, drop = FALSE]
  dplyr::bind_cols(ids, dplyr::bind_rows(c(list(prototype), results)))
}

#' @rdname dot-model_coef
#' @noRd
#' @export
.model_coef.Model <- function(model, n = Inf) {
  fits <- .model_fits(model)
  fits <- fits[seq_len(min(n, nrow(fits))), , drop = FALSE]
  result <- .bind_fit_results(
    fits,
    model$subject_id,
    lapply(fits$linear_model, .fit_coef),
    tibble::tibble(term = character(), estimate = double(), std.error = double(), df = double())
  )
  result[seq_len(min(n, nrow(result))), , drop = FALSE]
}

#' @rdname dot-model_coef
#' @noRd
#' @export
.model_coef.ModelFirth <- function(model, n = Inf) .model_coef.Model(model, n)

#' @rdname dot-model_coef
#' @noRd
#' @export
.model_coef.ModelLimma <- function(model, n = Inf) {
  # Moderation must use all subjects even when only a coefficient preview is requested.
  fit <- limma::eBayes(model$fit, trend = model$trend, robust = model$robust)
  rows <- seq_len(min(n, nrow(fit$coefficients)))
  result <- lapply(seq_len(min(n, ncol(fit$coefficients))), function(j) {
    dplyr::bind_cols(
      model$rowdata[rows, , drop = FALSE],
      tibble::tibble(
        term = colnames(fit$coefficients)[j],
        estimate = unname(fit$coefficients[rows, j]),
        std.error = unname(fit$stdev.unscaled[rows, j] * sqrt(fit$s2.post[rows])),
        df = unname(fit$df.total[rows])
      )
    )
  })
  result <- dplyr::bind_rows(result)
  result[seq_len(min(n, nrow(result))), , drop = FALSE]
}

#' @rdname dot-model_coef
#' @noRd
#' @export
.model_coef.NULL <- function(model, n = Inf) .model_not_available()

#' Coefficient covariance in the coefficient accessor's term order
#' @param fit a successful backend fit
#' @noRd
.fit_vcov <- function(fit) {
  terms <- .fit_coef(fit)$term
  as.matrix(stats::vcov(fit))[terms, terms, drop = FALSE]
}

#' Represent a named covariance matrix as a long table
#' @param covariance named coefficient covariance matrix
#' @noRd
.covariance_table <- function(covariance) {
  tibble::tibble(
    term = rep(rownames(covariance), times = ncol(covariance)),
    term2 = rep(colnames(covariance), each = nrow(covariance)),
    covariance = as.numeric(covariance)
  )
}

#' Coefficient covariance tables across subjects
#' @param model a model adapter
#' @noRd
.model_vcov <- function(model) UseMethod(".model_vcov")

#' @rdname dot-model_vcov
#' @noRd
#' @export
.model_vcov.Model <- function(model) {
  fits <- .model_fits(model)
  .bind_fit_results(
    fits,
    model$subject_id,
    lapply(fits$linear_model, function(fit) .covariance_table(.fit_vcov(fit))),
    tibble::tibble(term = character(), term2 = character(), covariance = double())
  )
}

#' @rdname dot-model_vcov
#' @noRd
#' @export
.model_vcov.ModelFirth <- function(model) .model_vcov.Model(model)

#' @rdname dot-model_vcov
#' @noRd
#' @export
.model_vcov.ModelLimma <- function(model) {
  fit <- limma::eBayes(model$fit, trend = model$trend, robust = model$robust)
  rows <- which(rowSums(!is.finite(fit$coefficients)) == 0L)
  result <- lapply(rows, function(i) {
    weights <- if (is.null(model$weights)) rep(1, ncol(model$E)) else model$weights[i, ]
    used <- is.finite(model$E[i, ]) & is.finite(weights) & weights > 0
    design <- fit$design[used, , drop = FALSE]
    covariance <- fit$s2.post[i] * solve(crossprod(design, weights[used] * design))
    .covariance_table(covariance)
  })
  .bind_fit_results(
    model$rowdata[rows, , drop = FALSE],
    model$subject_id,
    result,
    tibble::tibble(term = character(), term2 = character(), covariance = double())
  )
}

#' @rdname dot-model_vcov
#' @noRd
#' @export
.model_vcov.NULL <- function(model) .model_not_available()

#' Observation counts across model subjects
#' @param model a model adapter
#' @noRd
.model_nobs <- function(model) UseMethod(".model_nobs")

#' @rdname dot-model_nobs
#' @noRd
#' @export
.model_nobs.Model <- function(model) {
  fits <- .model_fits(model)
  dplyr::bind_cols(fits[, model$subject_id, drop = FALSE], nobs = vapply(fits$linear_model, .fit_nobs, numeric(1)))
}

#' @rdname dot-model_nobs
#' @noRd
#' @export
.model_nobs.ModelFirth <- function(model) .model_nobs.Model(model)

#' @rdname dot-model_nobs
#' @noRd
#' @export
.model_nobs.ModelLimma <- function(model) dplyr::bind_cols(model$rowdata, nobs = unname(rowSums(model$observed)))

#' @rdname dot-model_nobs
#' @noRd
#' @export
.model_nobs.NULL <- function(model) .model_not_available()

#' Model formula including its actual response and random effects
#' @param model a model adapter
#' @noRd
.model_formula <- function(model) UseMethod(".model_formula")

#' @rdname dot-model_formula
#' @noRd
#' @export
.model_formula.Model <- function(model) model$model_strategy$formula

#' @rdname dot-model_formula
#' @noRd
#' @export
.model_formula.ModelFirth <- function(model) model$models$strategy$formula

#' @rdname dot-model_formula
#' @noRd
#' @export
.model_formula.ModelLimma <- function(model) model$formula

#' @rdname dot-model_formula
#' @noRd
#' @export
.model_formula.NULL <- function(model) .model_not_available()

#' Model fit counts without evaluating any contrasts
#' @param model a model adapter
#' @noRd
.model_fit_status <- function(model) UseMethod(".model_fit_status")

#' Summarize the flags recorded for subject-wise fits
#' @param fits table of successful and failed subject fits
#' @noRd
.fit_status_table <- function(fits) {
  tibble::tibble(
    subjects = nrow(fits),
    fitted = sum(fits$has_model_fit),
    failed = sum(!fits$has_model_fit),
    singular = sum(fits$isSingular, na.rm = TRUE),
    imputed = sum(fits[["imputed"]], na.rm = TRUE),
    not_converged = sum(!fits[["converged"]], na.rm = TRUE)
  )
}

#' @rdname dot-model_fit_status
#' @noRd
#' @export
.model_fit_status.Model <- function(model) .fit_status_table(model$model_df)

#' @rdname dot-model_fit_status
#' @noRd
#' @export
.model_fit_status.ModelFirth <- function(model) .fit_status_table(model$models$model_df)

#' Identify matrix-model subjects with imputed responses
#' @param model a ModelLimma adapter
#' @noRd
.matrix_imputed_subjects <- function(model) {
  # limpa retains its imputed EList; model_name is only a user-facing label.
  if (!is.null(model$fit$EList)) {
    return(unname(rowSums(!model$observed) > 0))
  }
  model$rowdata[[model$subject_id[1]]] %in% model$imputed_proteins
}

#' @rdname dot-model_fit_status
#' @noRd
#' @export
.model_fit_status.ModelLimma <- function(model) {
  has_fit <- rowSums(is.na(model$fit$coefficients)) == 0
  tibble::tibble(
    subjects = nrow(model$rowdata),
    fitted = sum(has_fit),
    failed = sum(!has_fit),
    singular = NA_integer_,
    imputed = sum(.matrix_imputed_subjects(model)),
    not_converged = 0L
  )
}

#' @rdname dot-model_fit_status
#' @noRd
#' @export
.model_fit_status.NULL <- function(model) .model_not_available()

#' Fitted responses and residuals with observation identifiers
#' @param model a model adapter
#' @param lfqdata retained facade input providing column-role names
#' @noRd
.model_values <- function(model, lfqdata) UseMethod(".model_values")

#' Identifier columns to retain alongside observation-level results
#' @param model a model adapter
#' @param lfqdata retained facade input providing column-role names
#' @noRd
.value_keys <- function(model, lfqdata) {
  unique(c(model$subject_id, lfqdata$sample_name(), lfqdata$factor_keys(), lfqdata$hierarchy_keys()))
}

#' @rdname dot-model_values
#' @noRd
#' @export
.model_values.Model <- function(model, lfqdata) {
  fits <- .model_fits(model)
  keys <- setdiff(.value_keys(model, lfqdata), model$subject_id)
  if (nrow(fits) > 0) {
    keys <- intersect(keys, names(fits$data[[1]]))
  }
  results <- lapply(seq_len(nrow(fits)), function(i) {
    data <- fits$data[[i]]
    values <- .fit_values(fits$linear_model[[i]], data)
    dplyr::bind_cols(
      data[values$row, intersect(keys, names(data)), drop = FALSE],
      values[, c("response", "fitted", "residual")],
      imputed = isTRUE(fits[["imputed"]][i])
    )
  })
  source <- lfqdata$data_long()
  prototype <- dplyr::bind_cols(
    source[FALSE, intersect(keys, names(source)), drop = FALSE],
    tibble::tibble(response = double(), fitted = double(), residual = double(), imputed = logical())
  )
  .bind_fit_results(fits, model$subject_id, results, prototype)
}

#' @rdname dot-model_values
#' @noRd
#' @export
.model_values.ModelFirth <- function(model, lfqdata) .model_values.Model(model, lfqdata)

#' @rdname dot-model_values
#' @noRd
#' @export
.model_values.ModelLimma <- function(model, lfqdata) {
  fitted_matrix <- stats::fitted(model$fit)
  residual_matrix <- model$E - fitted_matrix
  observed <- model$observed & (rowSums(is.na(model$fit$coefficients)) == 0)
  cells <- which(observed, arr.ind = TRUE)
  keys <- intersect(setdiff(.value_keys(model, lfqdata), model$subject_id), names(model$annotation))
  dplyr::bind_cols(
    model$rowdata[cells[, 1], , drop = FALSE],
    model$annotation[cells[, 2], keys, drop = FALSE],
    tibble::tibble(
      response = model$E[cells],
      fitted = fitted_matrix[cells],
      residual = residual_matrix[cells],
      imputed = .matrix_imputed_subjects(model)[cells[, 1]]
    )
  )
}

#' @rdname dot-model_values
#' @noRd
#' @export
.model_values.NULL <- function(model, lfqdata) .model_not_available()

#' Fitted cell means and their uncertainty for one subject
#' @param fit a successful backend fit
#' @param contrast_fun the model strategy's contrast evaluator
#' @param level confidence level
#' @noRd
.fit_group_means <- function(fit, contrast_fun, level) UseMethod(".fit_group_means")

#' @rdname dot-fit_group_means
#' @noRd
#' @export
.fit_group_means.default <- function(fit, contrast_fun, level) {
  # Cell means align with observed factor combinations; marginal means do not.
  linear <- linfct_from_model(fit)$linfct_interactions
  estimates <- contrast_fun(fit, linear)
  margin <- stats::qt((1 + level) / 2, estimates$df) * estimates$std.error
  tibble::tibble(
    group = rownames(linear),
    estimate = estimates$estimate,
    conf.low = estimates$estimate - margin,
    conf.high = estimates$estimate + margin
  )
}

#' @rdname dot-fit_group_means
#' @noRd
#' @export
.fit_group_means.glm <- function(fit, contrast_fun, level) {
  result <- .fit_group_means.default(fit, contrast_fun, level)
  result[c("estimate", "conf.low", "conf.high")] <- lapply(
    result[c("estimate", "conf.low", "conf.high")],
    stats::plogis
  )
  result
}

#' @rdname dot-fit_group_means
#' @noRd
#' @export
.fit_group_means.logistf <- function(fit, contrast_fun, level) {
  data <- fit$model
  terms <- stats::terms(stats::formula(fit))
  factors <- intersect(attr(terms, "term.labels"), names(data))
  groups <- make_interaction_column(data, factors, sep = ":")$interaction
  rows <- which(!duplicated(groups))
  design <- stats::model.matrix(terms, data)[rows, , drop = FALSE]
  estimate <- as.numeric(design %*% stats::coef(fit))
  margin <- stats::qnorm((1 + level) / 2) * sqrt(rowSums((design %*% .fit_vcov(fit)) * design))
  tibble::tibble(
    group = as.character(groups[rows]),
    estimate = stats::plogis(estimate),
    conf.low = stats::plogis(estimate - margin),
    conf.high = stats::plogis(estimate + margin)
  )
}

#' Bind selected subjects' group means using the appropriate contrast strategy
#' @param model a per-subject model adapter
#' @param subjects identifiers in the first subject column
#' @param level confidence level
#' @param contrast_fun the model strategy's contrast evaluator
#' @noRd
.group_means_from_fits <- function(model, subjects, level, contrast_fun) {
  fits <- .model_fits(model)
  fits <- fits[fits[[model$subject_id[1]]] %in% subjects, , drop = FALSE]
  .bind_fit_results(
    fits,
    model$subject_id,
    lapply(fits$linear_model, .fit_group_means, contrast_fun = contrast_fun, level = level),
    tibble::tibble(group = character(), estimate = double(), conf.low = double(), conf.high = double())
  )
}

#' Group means and confidence limits for selected model subjects
#' @param model a model adapter
#' @param subjects identifiers in the first subject column
#' @param level confidence level
#' @noRd
.model_group_means <- function(model, subjects, level = 0.95) UseMethod(".model_group_means")

#' @rdname dot-model_group_means
#' @noRd
#' @export
.model_group_means.Model <- function(model, subjects, level = 0.95) {
  .group_means_from_fits(model, subjects, level, model$model_strategy$contrast_fun)
}

#' @rdname dot-model_group_means
#' @noRd
#' @export
.model_group_means.ModelFirth <- function(model, subjects, level = 0.95) {
  .group_means_from_fits(model, subjects, level, model$models$strategy$contrast_fun)
}

#' @rdname dot-model_group_means
#' @noRd
#' @export
.model_group_means.ModelLimma <- function(model, subjects, level = 0.95) {
  linear <- linfct_from_model(model$dummy_model)$linfct_interactions
  fit <- limma::contrasts.fit(model$fit, t(linear))
  fit <- limma::eBayes(fit, trend = model$trend, robust = model$robust)
  selected <- model$rowdata[[model$subject_id[1]]] %in% subjects & rowSums(!is.finite(model$fit$coefficients)) == 0L
  result <- lapply(seq_len(nrow(linear)), function(i) {
    table <- limma::topTable(fit, coef = i, number = Inf, sort.by = "none", confint = level)
    dplyr::bind_cols(
      model$rowdata[selected, , drop = FALSE],
      tibble::tibble(
        group = rep(rownames(linear)[i], sum(selected)),
        estimate = table$logFC[selected],
        conf.low = table$CI.L[selected],
        conf.high = table$CI.R[selected]
      )
    )
  })
  dplyr::bind_rows(result)
}

#' @rdname dot-model_group_means
#' @noRd
#' @export
.model_group_means.NULL <- function(model, subjects, level = 0.95) .model_not_available()
