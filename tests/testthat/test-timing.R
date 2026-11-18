# Fit time on noiseless data compared with equivalent noisy data (extended suite; timing is machine dependent).

#' @srrstats {RE7.1a} Extended test: lm and limma fit noiseless data (responses equal to the true group means) at
#'   least as fast as equivalent noisy data (median of 5 runs, ratio <= 1.2). rlm is excluded: on exact data its
#'   scale estimate is zero, its iterations do not converge and it is slower.
#' @noRd
NULL

test_that("lm and limma fit noiseless data at least as fast as noisy data (extended)", {
  skip_if_not_extended()
  sim <- sim_with_truth(Nprot = 500, N = 4, seed = 1)
  raw <- .with_seed(1, sim_lfq_data(Nprot = 500, N = 4))
  true_mean <- stats::setNames(raw$average_prot_abundance, paste(raw$proteinID, raw$idtype2, raw$sample, sep = "~"))
  data <- sim$lfqdata$data_long()
  data$abundance <- unname(true_mean[paste(data$protein_Id, data$sample, sep = "~")])
  expect_false(anyNA(data$abundance))
  noiseless <- LFQData$new(data, sim$lfqdata$get_config())
  fit_time <- function(method, lfqdata) {
    stats::median(replicate(
      5,
      system.time(
        muffle_warnings_matching(
          build_contrast_analysis(lfqdata, "~ group_", RECOVERY_CONTRASTS, method),
          "Zero sample variances detected|essentially perfect fit"
        )
      )[["elapsed"]]
    ))
  }
  for (method in c("lm", "limma")) {
    expect_lte(fit_time(method, noiseless), 1.2 * fit_time(method, sim$lfqdata))
  }
})
