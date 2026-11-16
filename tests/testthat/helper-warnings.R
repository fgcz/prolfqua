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

# plain limma on simulated data where some proteins miss a whole group: limma
# reports their NA coefficients; the *_impute methods refit exactly these
partial_na_warning <- "^Partial NA coefficients"

# proteins observed in one group only: their initial fit fails; the *_impute
# methods rescue exactly these, and plain methods report them in get_missing()
one_group_warning <- "contrasts can be applied only to factors with 2 or more levels"

# lmer on simulated proteins with a single peptide has no random-effect levels
single_peptide_warning <- "grouping factors must have > 1 sampled level"

# MASS::rlm on a few simulated proteins does not converge in its default 20 steps
rlm_convergence_warning <- "'rlm' failed to converge"
