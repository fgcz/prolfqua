#' @section Options:
#' \describe{
#'   \item{`prolfqua.verbose`}{`TRUE` (default) or `FALSE`. With `FALSE`,
#'     prolfqua prints no informative messages, such as columns added, rows
#'     removed or files written. Warnings and errors are not affected. The
#'     messages carry class `prolfqua_message`, so an application can route them
#'     to its own log with `withCallingHandlers(prolfqua_message = ...)`.}
#'   \item{`prolfqua.progress`}{Progress reporting while models and contrasts
#'     are computed per protein: `NULL` (default) shows a terminal progress
#'     bar, `"message"` a throttled message with elapsed time, and a
#'     `function(i, total, label)` is called with the progress.}
#' }
#' @keywords internal
"_PACKAGE"
