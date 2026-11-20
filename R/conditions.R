# Typed error conditions ----
#
# Internal helpers that raise `rlang` conditions carrying prolfqua-specific
# classes, so callers and tests can catch failures by class rather than by
# matching message text. All conditions inherit from `prolfqua_error`.

#' Abort because an argument fails its contract
#' @param arg name of the offending argument
#' @param must description of the requirement (e.g. "be a data frame")
#' @param not optional description of the actual value
#' @param parent optional parent condition to chain
#' @keywords internal
#' @noRd
abort_bad_argument <- function(arg, must, not = NULL, parent = NULL) {
  msg <- c("Invalid argument.", x = sprintf("`%s` must %s.", arg, must))
  if (!is.null(not)) {
    msg <- c(msg, i = sprintf("Actual value: %s.", not))
  }
  rlang::abort(
    msg,
    class = c("prolfqua_error_bad_argument", "prolfqua_error"),
    parent = parent
  )
}

#' Abort because required columns are missing from a data frame
#' @param cols character vector of missing column names
#' @param data_nm name to refer to the data frame by in the message
#' @keywords internal
#' @noRd
abort_missing_columns <- function(cols, data_nm = "data") {
  rlang::abort(
    c(
      "Required columns are missing.",
      x = sprintf("Missing from `%s`: %s.", data_nm, toString(cols))
    ),
    class = c("prolfqua_error_missing_column", "prolfqua_error_bad_argument", "prolfqua_error")
  )
}

#' Abort because the AnalysisConfiguration is invalid
#' @param msg one-line description of the problem
#' @keywords internal
#' @noRd
abort_invalid_config <- function(msg) {
  rlang::abort(
    c("Invalid configuration.", x = msg),
    class = c("prolfqua_error_invalid_configuration", "prolfqua_error")
  )
}

# Informative messages ----

#' Inform the user, unless silenced with `options(prolfqua.verbose = FALSE)`
#'
#' Signals an `rlang` message of class `prolfqua_message`, so applications can
#' route prolfqua's messages to their own log by class.
#' @param ... pasted together into the message, as in `message()`
#' @keywords internal
#' @noRd
inform_prolfqua <- function(...) {
  if (isTRUE(getOption("prolfqua.verbose", TRUE))) {
    rlang::inform(paste0(...), class = "prolfqua_message")
  }
  invisible(NULL)
}

#' Abort with a prolfqua error
#'
#' Like `stop()`: the arguments are pasted into the message. The condition has
#' class `prolfqua_error`, plus any classes given in `class`.
#' @param ... pasted together into the message
#' @param class additional condition classes, most specific first
#' @param call the call shown in the error message; defaults to the caller
#' @keywords internal
#' @noRd
abort_prolfqua <- function(..., class = NULL, call = rlang::caller_env()) {
  rlang::abort(paste0(...), class = c(class, "prolfqua_error"), call = call)
}

#' Match a case-sensitive choice exactly and retain prolfqua's error classes
#' @param arg character choice
#' @param values allowed choices
#' @param arg_nm argument name for the error
#' @noRd
#' @srrstats {G1.4a} All 107 previously undocumented top-level internal functions now have roxygen and @noRd blocks.
#' @srrstats {G2.3} Fixed choices are checked and builder help explicitly documents case sensitivity.
#' @srrstats {G2.3a} Exact arg_match is used here; legacy fixed-choice selectors use match.arg.
.arg_match_prolfqua <- function(arg, values, arg_nm) {
  if (!is.character(arg) || length(arg) != 1L || is.na(arg)) {
    abort_bad_argument(arg_nm, "be one character string")
  }
  tryCatch(
    rlang::arg_match(arg, values = values, error_arg = arg_nm),
    error = function(e) abort_bad_argument(arg_nm, paste0("be one of: ", toString(values)), parent = e)
  )
}

#' Parse a model formula with a classed argument error
#' @param modelstr one formula string or a formula
#' @noRd
#' @srrstatsTODO {G2.0} Formula and scalar lengths are checked here, and tidyMS_power_fdr.R checks all planning
#'   input lengths; limit: not all exported controls have been audited.
#' @srrstatsTODO {G2.6} Strategies accept formula objects and strings; limit: specialized vector classes are not fully
#'   tested.
.as_model_formula <- function(modelstr) {
  is_string <- is.character(modelstr) && length(modelstr) == 1 && !is.na(modelstr) && grepl("~", modelstr, fixed = TRUE)
  if (!(inherits(modelstr, "formula") || is_string)) {
    abort_bad_argument("modelstr", "be a formula or a single string containing `~`")
  }
  tryCatch(
    as.formula(modelstr),
    error = function(e) abort_bad_argument("modelstr", "be a valid model formula", parent = e)
  )
}

#' Require named contrast expressions
#' @param contrasts named character vector
#' @noRd
.assert_named_contrasts <- function(contrasts) {
  valid <- all(
    is.character(contrasts),
    length(contrasts) > 0,
    !anyNA(contrasts),
    !is.null(names(contrasts)),
    !anyNA(names(contrasts)),
    all(nzchar(names(contrasts)))
  )
  if (!valid) {
    abort_bad_argument("contrasts", "be a nonempty named character vector without missing values")
  }
}

#' Validate the data and configured response at input boundaries
#' @param data input data frame
#' @param config AnalysisConfiguration
#' @noRd
#' @srrstatsTODO {G2.1} LFQData/setup_analysis check data, configuration and response types; tidyMS_power_fdr.R
#'   rejects non-real numeric inputs and matrices; limit: the API audit is incomplete.
#' @srrstatsTODO {G2.4b} Numeric responses/counts are required and limma explicitly converts; no general conversion
#'   policy.
#' @srrstatsTODO {G2.12} Nonnumeric responses including lists are rejected; limit: other list columns lack a complete
#'   contract.
#' @srrstats {G5.8} Tests cover invalid types, designs, sparse inputs, empty inputs, non-finite responses, NA keys and
#'   all-NA or constant responses (test-input-validation.R, test-missing-nonfinite.R).
.assert_response_data <- function(data, config) {
  if (!is.data.frame(data)) {
    abort_bad_argument("data", "be a data frame")
  }
  if (!inherits(config, "AnalysisConfiguration")) {
    abort_invalid_config("`config` must be an AnalysisConfiguration.")
  }
  response <- config$get_response()
  if (length(response) != 1 || !response %in% names(data)) {
    abort_missing_columns(response %||% "response")
  }
  if (!is.numeric(data[[response]])) {
    abort_bad_argument("response", paste0("be numeric (column `", response, "`)"))
  }
}

#' Validate an LFQData before fitting or aggregation
#'
#' NA intensities are missing observations and pass. Zero rows, `Inf`, `-Inf` or `NaN` responses, NA in the
#' hierarchy, factor or sample key columns, and more than one row per sample and hierarchy key raise classed errors.
#' @param lfqdata LFQData to check
#' @noRd
#' @srrstats {EA2.1} Sample, hierarchy and isotope keys must identify each row; `setup_analysis()` and every model and
#'   aggregation entry (`.assert_model_input()`) check uniqueness and raise `prolfqua_error_duplicate_keys`.
.assert_model_input <- function(lfqdata) {
  data <- lfqdata$data_long()
  if (nrow(data) == 0) {
    abort_bad_argument("lfqdata", "have at least one row", not = "zero rows")
  }
  response <- lfqdata$response()
  values <- data[[response]]
  nonfinite <- !is.na(values) & !is.finite(values) | is.nan(values)
  if (any(nonfinite)) {
    subjects <- unique(do.call(paste, c(data[nonfinite, lfqdata$subject_id(), drop = FALSE], sep = "~")))
    abort_prolfqua(
      "Response `",
      response,
      "` has ",
      sum(nonfinite),
      " non-finite values (Inf, -Inf or NaN) in ",
      length(subjects),
      " subjects, e.g. ",
      toString(utils::head(subjects, 5)),
      ". Remove or replace them before modelling; NA marks a missing observation.",
      class = "prolfqua_error_nonfinite_response"
    )
  }
  keys <- unique(c(lfqdata$file_name(), lfqdata$sample_name(), lfqdata$factor_keys(), lfqdata$hierarchy_keys()))
  n_na <- vapply(keys, function(key) sum(is.na(data[[key]])), integer(1))
  if (any(n_na > 0)) {
    bad <- n_na[n_na > 0]
    abort_prolfqua(
      "Key columns contain NA: ",
      toString(paste0("`", names(bad), "` (", bad, ")")),
      ". Every row needs hierarchy, factor and sample values.",
      class = "prolfqua_error_missing_key"
    )
  }
  .assert_unique_keys(data, c(lfqdata$file_name(), lfqdata$hierarchy_keys(), lfqdata$isotope_label()))
  invisible(lfqdata)
}

#' Abort when the fixed-effect design is not of full rank across samples
#'
#' The check runs once per analysis on one row per sample, before any model is fitted. Random-effect terms
#' (`(1 | peptide_Id)`) are left out. Rank deficiency of single subjects caused by missing values is detected after
#' fitting instead (failed fit, listed by `get_missing()`).
#' @param lfqdata LFQData to fit
#' @param modelstr right-hand-side model formula string
#' @noRd
#' @srrstats {RE2.4, RE2.4a} Before fitting, every analysis checks that the fixed-effect design over the samples is of
#'   full rank and raises `prolfqua_error_collinear_design` naming the aliased coefficients (test-design-checks.R).
.assert_full_rank_design <- function(lfqdata, modelstr) {
  model_terms <- stats::terms(stats::as.formula(modelstr))
  labels <- attr(model_terms, "term.labels")
  labels <- labels[!grepl("|", labels, fixed = TRUE)]
  if (length(labels) == 0L) {
    return(invisible(lfqdata))
  }
  fixed <- stats::reformulate(labels, intercept = attr(model_terms, "intercept") == 1L)
  samples <- unique(lfqdata$data_long()[, unique(c(lfqdata$sample_name(), all.vars(fixed))), drop = FALSE])
  design <- stats::model.matrix(fixed, data = as.data.frame(samples))
  decomposition <- qr(design)
  if (decomposition$rank < ncol(design)) {
    aliased <- colnames(design)[decomposition$pivot[seq(decomposition$rank + 1L, ncol(design))]]
    abort_prolfqua(
      "The design `",
      paste(deparse(fixed), collapse = " "),
      "` is not of full rank across the ",
      nrow(samples),
      " samples: ",
      toString(paste0("`", aliased, "`")),
      " are linear combinations of the other columns (collinear or empty factor combinations). ",
      "Remove or merge the aliased factors.",
      class = "prolfqua_error_collinear_design"
    )
  }
  invisible(lfqdata)
}

#' Abort when key columns do not identify each row
#' @param data long-format data
#' @param keys sample, hierarchy and isotope key columns
#' @noRd
.assert_unique_keys <- function(data, keys) {
  duplicated_rows <- duplicated(data[, keys, drop = FALSE])
  if (!any(duplicated_rows)) {
    return(invisible(data))
  }
  ids <- unique(do.call(paste, c(data[duplicated_rows, keys, drop = FALSE], sep = "~")))
  abort_prolfqua(
    "More than one row for ",
    length(ids),
    " combinations of ",
    toString(paste0("`", keys, "`")),
    ", e.g. ",
    toString(utils::head(ids, 5)),
    ". Each hierarchy key must be observed at most once per sample; aggregate or deduplicate the rows first.",
    class = "prolfqua_error_duplicate_keys"
  )
}

#' Validate simulator size and missingness weight
#' @param Nprot positive whole number of proteins
#' @param weight_missing one number in [0, 1]
#' @noRd
#' @srrstatsTODO {G2.2} Scalar simulation arguments are enforced; limit: some thresholds and weights still need length
#'   checks.
.assert_simulation_args <- function(Nprot, weight_missing) {
  if (!.is_finite_number(Nprot) || !all(Nprot > 0, Nprot == floor(Nprot))) {
    abort_bad_argument("Nprot", "be one positive whole number")
  }
  if (!.is_finite_number(weight_missing) || !all(weight_missing >= 0, weight_missing <= 1)) {
    abort_bad_argument("weight_missing", "be one finite number in [0, 1]")
  }
}

#' Test whether a value is one finite real number
#' @param x value to test
#' @noRd
.is_finite_number <- function(x) {
  is.numeric(x) && length(x) == 1 && is.finite(x)
}


#' Merge user convergence controls over a backend's defaults
#'
#' @param control named list of controls given by the user
#' @param defaults named list of default controls (may be empty when the backend's own defaults apply)
#' @param allowed names the backend accepts
#' @return the merged named list
#' @noRd
.merge_control <- function(control, defaults, allowed = names(defaults)) {
  if (!is.list(control) || (length(control) > 0L && (is.null(names(control)) || any(names(control) == "")))) {
    abort_bad_argument("control", "be a named list")
  }
  unknown <- setdiff(names(control), allowed)
  if (length(unknown) > 0L) {
    abort_bad_argument(
      "control",
      paste0("only name ", if (length(allowed) > 0L) toString(allowed) else "nothing (no convergence controls)"),
      not = toString(unknown)
    )
  }
  utils::modifyList(defaults, control)
}
