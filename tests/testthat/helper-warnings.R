# Muffle only the warnings whose message matches `pattern` while evaluating
# `expr`; any other warning still reaches the test. Use it for warnings that
# are a known property of a test fixture, and name that property at the call.
muffle_warnings_matching <- function(expr, pattern) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(pattern, conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# loess of variance on the simulated peptide counts (1 to 8, mostly 1 to 3)
# has too few distinct counts near count 2 for a local quadratic fit
loess_ties_warning <- "pseudoinverse used|neighborhood radius|reciprocal condition number"
