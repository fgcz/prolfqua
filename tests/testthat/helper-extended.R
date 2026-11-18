# Extended tests: long-running recovery grids, nested-model recovery and timing comparisons.
# They run only when the environment variable PROLFQUA_EXTENDED_TESTS is "true"; see tests/README.md.

#' @srrstats {G5.10} Extended tests run in the same testthat suite and are switched on with
#'   PROLFQUA_EXTENDED_TESTS="true" (skip_if_not_extended()); CI sets the variable for commits whose message
#'   contains "run-extended".
#' @srrstats {G5.12} tests/README.md documents how to run the extended tests, their runtime, memory needs and that
#'   they download nothing and produce no artefacts for manual inspection.
#' @noRd
NULL

skip_if_not_extended <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("PROLFQUA_EXTENDED_TESTS"), "true"),
    "extended test; set PROLFQUA_EXTENDED_TESTS=true to run it"
  )
}
