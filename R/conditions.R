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
      x = sprintf("Missing from `%s`: %s.", data_nm, paste(cols, collapse = ", "))
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
