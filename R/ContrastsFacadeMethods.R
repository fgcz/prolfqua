#' Inspect a fitted prolfqua analysis
#'
#' Standard R methods for results returned by [build_contrast_analysis()].
#' Printing reports the model specification, fit status and a coefficient preview
#' without computing contrasts. Summaries evaluate the contrasts and report counts at the requested
#' FDR and absolute effect thresholds.
#'
#' @param x,object a fitted analysis returned by [build_contrast_analysis()].
#' @param ... arguments passed to the volcano or MA plotter by `plot()`;
#'   unused by the other methods.
#' @param FDR_threshold false discovery rate threshold for significant contrasts.
#' @param diff_threshold minimum absolute contrast effect for significance.
#' @param digits significant digits of the numbers printed by `print()` (default 3).
#' @return `print()` returns its input invisibly. `summary()` returns a list of
#'   class `summary.ContrastsFacadeBase`, with method, formula, fit status and a
#'   `per_contrast` table of estimated, missing, significant, up and down counts
#'   and median sigma and degrees of freedom when available.
#'   Matrix models report `NA` for singularity because it is not recorded by limma.
#'   `coef()` returns subject identifiers, term, estimate, standard error and df.
#'   `confint()` adds Wald t confidence limits (`conf.low`, `conf.high`).
#'   Limma uses moderated standard errors and degrees of freedom; Firth logistic
#'   fits use normal Wald intervals (infinite df). `formula()` returns the fitted
#'   model specification and `nobs()` returns a subject-keyed observation-count table.
#'   Counts exclude LOD-filled responses. Combined hurdle results have no single
#'   coefficient table, formula or observation count; inspect their intensity
#'   and detection facades through `analysis$contrast` instead.
#'   `fitted()` and `residuals()` return subject, sample, factor and available
#'   child identifiers, `response`, `fitted` or `residual`, and a subject-level
#'   `imputed` flag. Only observed responses from successful fits are included.
#'   Values use the fitting scale: transformed intensity or detection probability.
#'   Binomial responses are the fitted model's adjusted detection proportions;
#'   response residuals equal response minus fitted value. Hurdle results require
#'   selecting an intensity or detection component for these methods too.
#'   `plot()` returns a ggplot visibly. By default, observed responses and fitted
#'   cell means with confidence intervals are shown for six subjects selected by
#'   the smallest FDR in the first contrast. Mixed-model means use fixed effects
#'   and Satterthwaite uncertainty; detection means and limits are probabilities.
#'   ROPECA displays one panel per peptide of each selected protein. Hurdle
#'   analyses default to the volcano plot; request a component for fitted plots.
#'   `vcov()` returns subject identifiers, `term`, `term2` and `covariance`.
#'   It uses backend covariance (borrowed for LOD refits) or, for matrix models,
#'   the moderated variance times the inverse of the per-subject weighted design
#'   crossproduct. Failed matrix fits are omitted.
#'   `predict()` returns `newdata` with `prediction` on the fitting scale (probabilities for detection models;
#'   mixed models at population level, without random effects) and `imputed` (the subject's fit is a
#'   limit-of-detection refit). Hurdle results predict per component.
#'   The fit status counts subjects that did not converge (`not_converged`); per subject, `converged` and
#'   `iterations` are columns of the model table (`analysis$model$model_df`).
#' @name prolfqua_analysis_methods
#' @family modelling
#' @examples
#' sim <- sim_lfq_data_protein_config(n_proteins = 10, weight_missing = 0, seed = 42)
#' lfq <- LFQData$new(sim$data, sim$config)
#' analysis <- build_contrast_analysis(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"))
#' print(analysis)
#' summary(analysis)
#' head(stats::coef(analysis))
#' head(stats::vcov(analysis))
#' head(stats::confint(analysis))
#' stats::formula(analysis)
#' head(stats::nobs(analysis))
#' head(stats::fitted(analysis))
#' head(stats::residuals(analysis))
#' plot(analysis)
#' plot(analysis, type = "volcano")
NULL

#' Formula information for a single model or both hurdle components
#' @param object a fitted analysis facade
#' @noRd
.facade_formulas <- function(object) {
  if (is.null(object$model)) {
    return(list(
      intensity = .model_formula(object$contrast$intensity$model),
      detection = .model_formula(object$contrast$detection$model)
    ))
  }
  .model_formula(object$model)
}

#' Fit counts for a single model or both hurdle components
#' @param object a fitted analysis facade
#' @noRd
.facade_fit_status <- function(object) {
  if (is.null(object$model)) {
    return(dplyr::bind_rows(
      intensity = .model_fit_status(object$contrast$intensity$model),
      detection = .model_fit_status(object$contrast$detection$model),
      .id = "component"
    ))
  }
  .model_fit_status(object$model)
}

#' Write fitted formula information without printing its environment
#' @param formulas a formula or a named list of hurdle component formulas
#' @noRd
.print_analysis_formulas <- function(formulas) {
  if (inherits(formulas, "formula")) {
    cat("Formula:", paste(deparse(formulas), collapse = " "), "\n")
  } else {
    for (name in names(formulas)) {
      cat(name, "formula:", paste(deparse(formulas[[name]]), collapse = " "), "\n")
    }
  }
}

#' @rdname prolfqua_analysis_methods
#' @exportS3Method base::print
#' @srrstats {RE4.17} Printing shows the model specification, fit counts and coefficients without evaluating contrasts.
#' @srrstats {EA4.1, EA5.2} `print()` of analyses, their summaries and LFQData rounds printed numbers to `digits`
#'   significant digits (default 3); `contrast_summary_table(digits = )` does the same for report tables.
print.ContrastsFacadeBase <- function(x, digits = 3, ...) {
  cat("prolfqua analysis:", x$facade_name, "(", class(x)[1], ")\n")
  .print_analysis_formulas(.facade_formulas(x))
  cat("Contrasts:", paste(x$.contrast_names, collapse = ", "), "\n")
  input <- x$.lfqdata
  cat("Input subjects:", nrow(unique(input$data_long()[, input$subject_id(), drop = FALSE])), "\n")
  print(.facade_fit_status(x))
  cat(.provenance_line(x$.recorded_provenance()), "\n")
  if (is.null(x$model)) {
    cat("Intensity coefficient preview:\n")
    print(.signif_columns(.model_coef(x$contrast$intensity$model, n = 6L), digits))
    cat("Detection coefficient preview:\n")
    print(.signif_columns(.model_coef(x$contrast$detection$model, n = 6L), digits))
  } else {
    cat("Coefficient preview:\n")
    print(.signif_columns(.model_coef(x$model, n = 6L), digits))
  }
  invisible(x)
}

#' Round the double columns of a table to significant digits for printing
#' @param data data frame
#' @param digits significant digits
#' @noRd
.signif_columns <- function(data, digits) {
  doubles <- vapply(data, is.double, logical(1))
  data[doubles] <- lapply(data[doubles], signif, digits = digits)
  data
}

#' Median of an optional contrast statistic
#' @param data contrast result rows
#' @param column optional numeric column name
#' @noRd
.contrast_median <- function(data, column) {
  if (!column %in% names(data)) {
    return(NA_real_)
  }
  stats::median(data[[column]], na.rm = TRUE)
}

#' @rdname prolfqua_analysis_methods
#' @exportS3Method base::summary
#' @srrstats {RE4.18} Summaries compute contrast statistics separately from printing the fitted object.
summary.ContrastsFacadeBase <- function(object, FDR_threshold = 0.05, diff_threshold = 1, ...) {
  contrasts <- object$get_contrasts()
  missing <- object$get_missing()
  significant <- object$filter_significant(FDR_threshold = FDR_threshold, diff_threshold = diff_threshold)
  config <- object$get_config()
  per_contrast <- lapply(object$.contrast_names, function(name) {
    rows <- contrasts[contrasts[[config$contrast_col]] == name, , drop = FALSE]
    selected <- significant[significant[[config$contrast_col]] == name, , drop = FALSE]
    tibble::tibble(
      contrast = name,
      estimated = nrow(rows),
      missing = sum(missing[[config$contrast_col]] == name),
      significant = nrow(selected),
      up = sum(selected[[config$effect_col]] > 0),
      down = sum(selected[[config$effect_col]] < 0),
      median_sigma = .contrast_median(rows, "sigma"),
      median_df = .contrast_median(rows, "df")
    )
  })
  structure(
    list(
      method = object$facade_name,
      formula = .facade_formulas(object),
      fit_status = .facade_fit_status(object),
      per_contrast = dplyr::bind_rows(per_contrast),
      provenance = .provenance_line(object$provenance())
    ),
    class = "summary.ContrastsFacadeBase"
  )
}

#' @rdname prolfqua_analysis_methods
#' @exportS3Method base::print
print.summary.ContrastsFacadeBase <- function(x, digits = 3, ...) {
  cat("prolfqua analysis summary:", x$method, "\n")
  .print_analysis_formulas(x$formula)
  print(x$fit_status)
  print(.signif_columns(x$per_contrast, digits))
  cat(x$provenance, "\n")
  invisible(x)
}

#' @rdname prolfqua_analysis_methods
#' @exportS3Method stats::coef
coef.ContrastsFacadeBase <- function(object, ...) .model_coef(object$model)

#' @rdname prolfqua_analysis_methods
#' @exportS3Method stats::vcov
#' @srrstats {RE4.6} Subject-keyed covariance tables expose backend, borrowed and moderated matrix covariances.
vcov.ContrastsFacadeBase <- function(object, ...) .model_vcov(object$model)

#' Resolve requested coefficient names or indices
#' @param parm requested coefficient names or positive indices
#' @param terms available coefficient names
#' @noRd
.coefficient_terms <- function(parm, terms) {
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(!is.finite(parm) | parm < 1 | parm > length(terms) | parm != floor(parm))) {
      abort_bad_argument("parm", "contain valid positive coefficient indices")
    }
    parm <- terms[parm]
  }
  if (!is.character(parm) || anyNA(parm) || !all(parm %in% terms)) {
    abort_bad_argument("parm", "contain coefficient names from coef(object)")
  }
  parm
}

#' @rdname prolfqua_analysis_methods
#' @param parm coefficient names or positive coefficient indices; omitted for all terms.
#' @param level confidence level, strictly between zero and one.
#' @exportS3Method stats::confint
#' @srrstats {RE4.3} Coefficient confidence intervals use Wald t limits, moderated limma uncertainty and normal Firth
#'   limits.
confint.ContrastsFacadeBase <- function(object, parm, level = 0.95, ...) {
  if (!.is_finite_number(level) || level <= 0 || level >= 1) {
    abort_bad_argument("level", "be one finite number strictly between zero and one")
  }
  coefficients <- .model_coef(object$model)
  if (!missing(parm)) {
    terms <- unique(coefficients$term)
    parm <- .coefficient_terms(parm, terms)
    coefficients <- coefficients[coefficients$term %in% parm, , drop = FALSE]
  }
  margin <- stats::qt(1 - (1 - level) / 2, df = coefficients$df) * coefficients$std.error
  coefficients$conf.low <- coefficients$estimate - margin
  coefficients$conf.high <- coefficients$estimate + margin
  coefficients
}

#' @rdname prolfqua_analysis_methods
#' @exportS3Method stats::formula
#' @srrstats {RE4.4} The accessor returns the fitted specification, including nested and detection-model terms.
formula.ContrastsFacadeBase <- function(x, ...) .model_formula(x$model)

#' @rdname prolfqua_analysis_methods
#' @exportS3Method stats::nobs
#' @srrstats {RE4.5} Subject-keyed observation counts exclude imputed responses and unavailable matrix cells.
nobs.ContrastsFacadeBase <- function(object, ...) .model_nobs(object$model)

#' @rdname prolfqua_analysis_methods
#' @exportS3Method stats::fitted
#' @srrstats {RE4.9} Fitted responses retain subject, sample, factor and child identifiers on the documented fitting
#'   scale.
fitted.ContrastsFacadeBase <- function(object, ...) {
  values <- .model_values(object$model, object$.lfqdata)
  values[, setdiff(names(values), "residual"), drop = FALSE]
}

#' @rdname prolfqua_analysis_methods
#' @param newdata for `predict()`: data frame with the subject identifier and predictor columns; `NULL` (default)
#'   predicts every subject in every sample of the analysis input, including samples with missing responses.
#' @exportS3Method stats::predict
#' @srrstats {RE4.16} `predict()` accepts any combination of subject and group, also groups absent from a subject's
#'   fit: methods with a limit-of-detection refit (`_impute`) predict them from the refit, others give NA. Limit: a
#'   group absent from the whole experiment has no coefficient. Tested in test-predict.R.
predict.ContrastsFacadeBase <- function(object, newdata = NULL, ...) {
  if (is.null(object$model)) {
    .model_not_available()
  }
  subject_id <- object$model$subject_id
  if (is.null(newdata)) {
    input <- object$.lfqdata
    subjects <- dplyr::distinct(input$data_long()[, subject_id, drop = FALSE])
    newdata <- tidyr::crossing(subjects, .sample_template(input))
  }
  if (!is.data.frame(newdata)) {
    abort_bad_argument("newdata", "be a data frame or NULL")
  }
  missing <- setdiff(subject_id, names(newdata))
  if (length(missing) > 0L) {
    abort_missing_columns(missing)
  }
  .model_predict(object$model, newdata)
}

#' @rdname prolfqua_analysis_methods
#' @exportS3Method stats::residuals
#' @srrstats {RE4.10} Observed response residuals and their fitting scale are documented and tested.
residuals.ContrastsFacadeBase <- function(object, ...) {
  values <- .model_values(object$model, object$.lfqdata)
  values[, setdiff(names(values), "fitted"), drop = FALSE]
}

#' @rdname prolfqua_analysis_methods
#' @param type plot type: fitted responses (default), FDR volcano or MA.
#'   Hurdle analyses default to volcano because they have two models.
#' @param subject character vector of identifiers from the first model subject
#'   column for fitted plots; `NULL` selects up to six by FDR in the first contrast.
#' @exportS3Method graphics::plot
#' @srrstats {RE6.0} The default plot shows fitted responses, or an FDR volcano for combined hurdle analyses.
#' @srrstats {RE6.1} S3 plot dispatch connects analysis results to the existing contrast plotter.
#' @srrstats {RE6.2} Default fitted plots show observations and group means with configurable confidence limits.
plot.ContrastsFacadeBase <- function(x, type = c("fitted", "volcano", "ma"), subject = NULL, level = 0.95, ...) {
  if (missing(type)) {
    type <- if (is.null(x$model)) "volcano" else "fitted"
  }
  type <- .arg_match_prolfqua(type, c("fitted", "volcano", "ma"), "type")
  if (type == "fitted") {
    return(.plot_fitted_analysis(x, subject, level))
  }
  plotter <- x$get_Plotter()
  if (type == "volcano") {
    plotter$volcano(...)$FDR
  } else {
    plotter$ma_plot(...)
  }
}

#' Resolve explicitly requested subjects or the six smallest first-contrast FDRs
#' @param x a fitted analysis facade
#' @param subject requested identifiers or NULL for selection by FDR
#' @noRd
.fitted_plot_subjects <- function(x, subject) {
  available <- .model_nobs(x$model)
  key <- x$model$subject_id[1]
  if (is.null(subject)) {
    contrasts <- x$get_contrasts()
    config <- x$get_config()
    contrasts <- contrasts[contrasts[[config$contrast_col]] == x$.contrast_names[1], , drop = FALSE]
    ordered <- contrasts[[key]][order(contrasts[[config$fdr_col]])]
    subject <- head(unique(as.character(ordered[ordered %in% available[[key]]])), 6L)
  }
  if (!is.character(subject) || length(subject) == 0L || anyNA(subject) || !all(subject %in% available[[key]])) {
    abort_bad_argument("subject", "contain one or more identifiers from the first model subject column")
  }
  unique(subject)
}

#' Plot observed responses and estimated cell means on the fitted response scale
#' @param x a fitted analysis facade
#' @param subject selected first-column subject identifiers or NULL
#' @param level confidence level
#' @noRd
.plot_fitted_analysis <- function(x, subject, level) {
  if (!.is_finite_number(level) || level <= 0 || level >= 1) {
    abort_bad_argument("level", "be one finite number strictly between zero and one")
  }
  subject <- .fitted_plot_subjects(x, subject)
  keys <- x$model$subject_id
  observed <- .model_values(x$model, x$.lfqdata)
  observed <- observed[observed[[keys[1]]] %in% subject, , drop = FALSE]
  factors <- intersect(attr(stats::terms(.model_formula(x$model)), "term.labels"), x$.lfqdata$factor_keys())
  observed$group <- as.character(make_interaction_column(observed, factors, sep = ":")$interaction)
  means <- .model_group_means(x$model, subject, level)
  observed$subject <- do.call(paste, c(observed[keys], sep = " / "))
  means$subject <- do.call(paste, c(means[keys], sep = " / "))
  ggplot2::ggplot(observed, ggplot2::aes(x = .data$group, y = .data$response)) +
    theme_prolfqua(grid = TRUE) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1), alpha = 0.5) +
    ggplot2::geom_pointrange(
      data = means,
      mapping = ggplot2::aes(x = .data$group, y = .data$estimate, ymin = .data$conf.low, ymax = .data$conf.high),
      inherit.aes = FALSE,
      colour = "#0072B2"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$subject), scales = "free_y") +
    ggplot2::labs(x = "Group", y = "Modelled response")
}


#' One-line provenance summary for print() and summary()
#' @param entries provenance entries
#' @noRd
.provenance_line <- function(entries) {
  paste0(
    "Provenance: ",
    sum(entries$failure),
    " failure and ",
    sum(!entries$failure),
    " other entries; see $provenance()"
  )
}
