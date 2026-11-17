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
