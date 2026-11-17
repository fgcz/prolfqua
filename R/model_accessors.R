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
#' @noRd
.model_coef <- function(model) UseMethod(".model_coef")

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
.model_coef.Model <- function(model) {
  fits <- .model_fits(model)
  .bind_fit_results(
    fits,
    model$subject_id,
    lapply(fits$linear_model, .fit_coef),
    tibble::tibble(term = character(), estimate = double(), std.error = double(), df = double())
  )
}

#' @rdname dot-model_coef
#' @noRd
#' @export
.model_coef.ModelFirth <- function(model) .model_coef.Model(model)

#' @rdname dot-model_coef
#' @noRd
#' @export
.model_coef.ModelLimma <- function(model) {
  fit <- limma::eBayes(model$fit, trend = model$trend, robust = model$robust)
  result <- lapply(seq_len(ncol(fit$coefficients)), function(j) {
    dplyr::bind_cols(
      model$rowdata,
      tibble::tibble(
        term = colnames(fit$coefficients)[j],
        estimate = unname(fit$coefficients[, j]),
        std.error = unname(fit$stdev.unscaled[, j] * sqrt(fit$s2.post)),
        df = unname(fit$df.total)
      )
    )
  })
  dplyr::bind_rows(result)
}

#' @rdname dot-model_coef
#' @noRd
#' @export
.model_coef.NULL <- function(model) .model_not_available()

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
    imputed = sum(fits[["imputed"]], na.rm = TRUE)
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
  fitted <- rowSums(is.na(model$fit$coefficients)) == 0
  tibble::tibble(
    subjects = nrow(model$rowdata),
    fitted = sum(fitted),
    failed = sum(!fitted),
    singular = sum(!fitted),
    imputed = sum(.matrix_imputed_subjects(model))
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
