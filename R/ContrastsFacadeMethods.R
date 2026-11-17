#' Inspect a fitted prolfqua analysis
#'
#' Standard R methods for results returned by [build_contrast_analysis()].
#' Printing reports the model specification and fit status without computing
#' contrasts. Summaries evaluate the contrasts and report counts at the requested
#' FDR and absolute effect thresholds.
#'
#' @param x,object a fitted analysis returned by [build_contrast_analysis()].
#' @param ... additional arguments; unused by print and summary.
#' @param FDR_threshold false discovery rate threshold for significant contrasts.
#' @param diff_threshold minimum absolute contrast effect for significance.
#' @return `print()` returns its input invisibly. `summary()` returns a list of
#'   class `summary.ContrastsFacadeBase`, with method, formula, fit status and a
#'   `per_contrast` table of estimated, missing, significant, up and down counts
#'   and median sigma and degrees of freedom when available.
#'   `coef()` returns subject identifiers, term, estimate, standard error and df.
#'   `confint()` adds Wald t confidence limits (`conf.low`, `conf.high`).
#'   Limma uses moderated standard errors and degrees of freedom; Firth logistic
#'   fits use normal Wald intervals (infinite df). `formula()` returns the fitted
#'   model specification and `nobs()` returns a subject-keyed observation-count table.
#'   Counts exclude LOD-filled responses. Combined hurdle results have no single
#'   coefficient table, formula or observation count; inspect their intensity
#'   and detection facades through `analysis$contrast` instead.
#' @name prolfqua_analysis_methods
#' @family modelling
#' @examples
#' sim <- sim_lfq_data_protein_config(Nprot = 10, weight_missing = 0, seed = 42)
#' lfq <- LFQData$new(sim$data, sim$config)
#' analysis <- build_contrast_analysis(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"))
#' print(analysis)
#' summary(analysis)
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
print.ContrastsFacadeBase <- function(x, ...) {
  cat("prolfqua analysis:", x$facade_name, "(", class(x)[1], ")\n")
  .print_analysis_formulas(.facade_formulas(x))
  cat("Contrasts:", paste(x$.contrast_names, collapse = ", "), "\n")
  input <- x$.lfqdata
  cat("Input subjects:", nrow(unique(input$data_long()[, input$subject_id(), drop = FALSE])), "\n")
  print(.facade_fit_status(x))
  invisible(x)
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
      per_contrast = dplyr::bind_rows(per_contrast)
    ),
    class = "summary.ContrastsFacadeBase"
  )
}

#' @rdname prolfqua_analysis_methods
#' @exportS3Method base::print
print.summary.ContrastsFacadeBase <- function(x, ...) {
  cat("prolfqua analysis summary:", x$method, "\n")
  .print_analysis_formulas(x$formula)
  print(x$fit_status)
  print(x$per_contrast)
  invisible(x)
}

#' @rdname prolfqua_analysis_methods
#' @exportS3Method stats::coef
coef.ContrastsFacadeBase <- function(object, ...) .model_coef(object$model)

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
formula.ContrastsFacadeBase <- function(x, ...) .model_formula(x$model)

#' @rdname prolfqua_analysis_methods
#' @exportS3Method stats::nobs
nobs.ContrastsFacadeBase <- function(object, ...) .model_nobs(object$model)
