# Parameter recovery on simulated data with known group differences.

#' @srrstats {G5.6} Contrast estimates of all aggregated intensity methods recover simulated true group differences
#'   (bias, slope of estimate on truth, RMSE within tolerance); nested methods recover the true peptide-level
#'   differences in the extended suite.
#' @srrstats {G5.6b} Recovery is repeated over seeds 1-2 in the standard suite and 1-50 in the extended suite, with
#'   tolerances per seed and tighter ones on the average over seeds.
#' @srrstats {G5.7} Estimation error falls as the number of samples per group grows from 3 to 6 to 12 (and to 24 in
#'   the extended suite).
#' @noRd
NULL

AGGREGATED_RECOVERY_METHODS <- c(
  "lm",
  "lm_impute",
  "rlm",
  "rfit",
  "limma",
  "limma_impute",
  "limma_voom",
  "deqms",
  "deqms_voom"
)
RECOVERY_WARNINGS <- paste(rlm_convergence_warning, loess_ties_warning, single_peptide_warning, sep = "|")

fit_contrasts <- function(method, lfqdata, modelstr = "~ group_") {
  muffle_warnings_matching(
    build_contrast_analysis(lfqdata, modelstr, RECOVERY_CONTRASTS, method)$get_contrasts(),
    RECOVERY_WARNINGS
  )
}

# Per seed: |bias| < 0.25, |slope - 1| < 0.3, RMSE < max_rmse. Over seeds: |mean bias| < 0.08, |mean slope - 1| < 0.1.
expect_recovery <- function(method, seeds, N = 4) {
  stats <- do.call(
    rbind,
    lapply(seeds, function(seed) {
      sim <- sim_with_truth(n_proteins = 200, n_per_group = N, seed = seed, counts = grepl("deqms", method))
      recovery_stats(with_truth(fit_contrasts(method, sim$lfqdata), sim$truth))
    })
  )
  expect_true(all(abs(stats$bias) < 0.25), label = paste(method, "bias per seed"))
  expect_true(all(abs(stats$slope - 1) < 0.3), label = paste(method, "slope per seed"))
  expect_true(all(stats$rmse < 0.9), label = paste(method, "RMSE per seed"))
  expect_lt(abs(mean(stats$bias)), 0.08)
  expect_lt(abs(mean(stats$slope) - 1), 0.1)
  invisible(stats)
}

skip_if_backend_missing <- function(method) {
  if (method == "rfit") skip_if_not_installed("Rfit")
}

for (method in AGGREGATED_RECOVERY_METHODS) {
  test_that(paste(method, "recovers simulated fold changes over 2 seeds"), {
    skip_if_backend_missing(method)
    expect_recovery(method, seeds = 1:2)
  })
  test_that(paste(method, "recovers simulated fold changes over 50 seeds (extended)"), {
    skip_if_not_extended()
    skip_if_backend_missing(method)
    expect_recovery(method, seeds = 1:50)
  })
}

rmse_by_group_size <- function(method, sizes, seeds) {
  vapply(
    sizes,
    function(N) {
      mean(vapply(
        seeds,
        function(seed) {
          sim <- sim_with_truth(n_proteins = 200, n_per_group = N, seed = seed)
          mean(recovery_stats(with_truth(fit_contrasts(method, sim$lfqdata), sim$truth))$rmse)
        },
        numeric(1)
      ))
    },
    numeric(1)
  )
}

test_that("estimation error falls as samples per group grow", {
  for (method in c("lm", "limma")) {
    rmse <- rmse_by_group_size(method, sizes = c(3, 6, 12), seeds = 1:3)
    expect_true(all(diff(rmse) < 0), label = paste(method, "RMSE decreases"))
    # The standard error of a difference of two means scales with 1 / sqrt(N): expected ratio 0.5.
    expect_lt(rmse[3] / rmse[1], 0.6)
  }
})

test_that("estimation error falls with group size for all aggregated methods (extended)", {
  skip_if_not_extended()
  for (method in AGGREGATED_RECOVERY_METHODS) {
    skip_if_backend_missing(method)
    rmse <- rmse_by_group_size(method, sizes = c(3, 6, 12, 24), seeds = 1:10)
    expect_true(all(diff(rmse) < 0), label = paste(method, "RMSE decreases"))
  }
})

# Nested methods estimate the protein effect from peptide data: lmer and limpa a mean over peptides,
# ROPECA the median of the peptide estimates.
NESTED_RECOVERY <- list(
  lmer_nested = list(truth = "peptide_mean", modelstr = "~ group_ + (1 | peptide_Id) + (1 | sampleName)", slope = 0.1),
  limpa_nested = list(truth = "peptide_mean", modelstr = "~ group_", slope = 0.25),
  ropeca_nested = list(truth = "peptide_median", modelstr = "~ group_", slope = 0.1),
  ropecaV2_nested = list(truth = "peptide_median", modelstr = "~ group_", slope = 0.1)
)
NESTED_PACKAGES <- list(lmer_nested = c("lme4", "lmerTest"), limpa_nested = "limpa")

test_that("nested methods recover the simulated peptide-level differences (extended)", {
  skip_if_not_extended()
  for (method in names(NESTED_RECOVERY)) {
    for (package in NESTED_PACKAGES[[method]]) {
      skip_if_not_installed(package)
    }
    spec <- NESTED_RECOVERY[[method]]
    stats <- do.call(
      rbind,
      lapply(1:20, function(seed) {
        sim <- sim_with_truth(n_proteins = 100, n_per_group = 4, seed = seed, peptide = TRUE)
        recovery_stats(with_truth(fit_contrasts(method, sim$lfqdata, spec$modelstr), sim$truth, spec$truth))
      })
    )
    expect_lt(abs(mean(stats$bias)), 0.1)
    expect_lt(abs(mean(stats$slope) - 1), spec$slope)
  }
})
