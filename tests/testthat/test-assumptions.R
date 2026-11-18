# Implications of violated input-data assumptions, as documented in ?prolfqua (Glossary,
# Linear-model assumptions) and ?LFQData (Missing Data Assumptions).

#' @srrstats {RE1.4} Assumptions on input data are documented in ?prolfqua (normal, homoscedastic errors) and
#'   ?LFQData (missing at random); tests show the documented implications of violating them: outliers pull lm
#'   estimates while rlm and rfit resist them, and intensity-dependent missingness shrinks fold changes of
#'   down-regulated proteins toward zero, also after LOD imputation.
#' @noRd
NULL

# Removing low intensities leaves some proteins with one group; their fits fail as expected.
ASSUMPTION_WARNINGS <- paste(rlm_convergence_warning, partial_na_warning, one_group_warning, sep = "|")

assumption_contrasts <- function(method, lfqdata) {
  muffle_warnings_matching(
    build_contrast_analysis(lfqdata, "~ group_", RECOVERY_CONTRASTS, method)$get_contrasts(),
    ASSUMPTION_WARNINGS
  )
}

test_that("outliers pull lm estimates while rlm and rfit resist them", {
  skip_if_not_installed("Rfit")
  for (seed in 1:3) {
    sim <- sim_with_truth(Nprot = 200, N = 4, seed = seed)
    data <- sim$lfqdata$data_long()
    # Shift 5% of the observations by 6 units, six residual standard deviations.
    outliers <- .with_seed(100 + seed, sample(nrow(data), round(0.05 * nrow(data))))
    data$abundance[outliers] <- data$abundance[outliers] + 6
    contaminated <- LFQData$new(data, sim$lfqdata$get_config())
    rmse <- vapply(
      c("lm", "rlm", "rfit"),
      function(method) {
        joined <- with_truth(assumption_contrasts(method, contaminated), sim$truth)
        sqrt(mean((joined$diff - joined$true_diff)^2))
      },
      numeric(1)
    )
    expect_lt(rmse[["rlm"]], 0.9 * rmse[["lm"]])
    expect_lt(rmse[["rfit"]], 0.9 * rmse[["lm"]])
  }
})

test_that("missing low intensities shrink fold changes of down-regulated proteins toward zero", {
  attenuation <- numeric()
  imputed_bias <- numeric()
  for (seed in 1:5) {
    sim <- sim_with_truth(Nprot = 200, N = 4, seed = seed)
    data <- sim$lfqdata$data_long()
    # Missing not at random: the lowest 15% of all intensities are not observed.
    observed <- data[data$abundance > stats::quantile(data$abundance, 0.15), ]
    mnar <- LFQData$new(observed, sim$lfqdata$get_config())
    complete <- with_truth(assumption_contrasts("lm", sim$lfqdata), sim$truth)
    complete_case <- with_truth(assumption_contrasts("lm", mnar), sim$truth)
    matched <- merge(complete, complete_case, by = c("protein_Id", "contrast", "true_diff"), suffixes = c("", "_mnar"))
    down <- matched$true_diff < -0.5
    # Same proteins with and without the missing values: estimates move toward zero.
    attenuation <- c(attenuation, mean(matched$diff_mnar[down] - matched$diff[down]))
    imputed <- with_truth(assumption_contrasts("lm_impute", mnar), sim$truth)
    rescued <- !paste(imputed$protein_Id, imputed$contrast) %in% paste(complete_case$protein_Id, complete_case$contrast)
    rescued <- rescued & imputed$true_diff < -0.5
    # Proteins estimable only after imputing at the LOD: the LOD caps how low they can go.
    imputed_bias <- c(imputed_bias, mean(imputed$diff[rescued] - imputed$true_diff[rescued]))
  }
  expect_true(all(attenuation > 0))
  expect_gt(mean(attenuation), 0.05)
  expect_true(all(imputed_bias > 0.5))
})
