# Tests for ropecaV2_nested: ROPECA protein p-values with an exact null for the
# median peptide z-score and the correlation of peptides within a protein.

# Median of n peptide z-scores per protein; the peptides of a protein share a
# protein effect, so their z-scores have correlation rho.
simulate_median_z <- function(n_prot, n, rho, seed) {
  withr::local_seed(seed)
  shared <- stats::rnorm(n_prot)
  own <- matrix(stats::rnorm(n_prot * n), n_prot)
  z <- sqrt(rho) * shared + sqrt(1 - rho) * own
  apply(z, 1, stats::median)
}

# Peptide abundances, 4 vs 4 samples; peptides of a protein share a sample effect
# with variance rho (total variance 1), some values missing at random.
simulate_correlated_peptides <- function(n_prot, n_pep, rho, missing = 0, seed = 1) {
  withr::local_seed(seed)
  samples <- paste0(rep(c("A", "B"), each = 4), "_", 1:4)
  data <- expand.grid(
    sampleName = samples,
    peptide = seq_len(n_pep),
    protein = seq_len(n_prot),
    stringsAsFactors = FALSE
  )
  sample_effect <- matrix(stats::rnorm(n_prot * length(samples), sd = sqrt(rho)), n_prot)
  data$abundance <- 20 +
    sample_effect[cbind(data$protein, match(data$sampleName, samples))] +
    stats::rnorm(nrow(data), sd = sqrt(1 - rho))
  data$abundance[stats::runif(nrow(data)) < missing] <- NA
  data.frame(
    protein_Id = paste0("P", data$protein),
    peptide_Id = paste0("P", data$protein, "_", data$peptide),
    sampleName = data$sampleName,
    group_ = substr(data$sampleName, 1, 1),
    abundance = data$abundance
  )
}

test_that(".median_normal_sf matches simulated medians of normals", {
  withr::local_seed(1)
  for (n in 1:6) {
    medians <- apply(matrix(stats::rnorm(4e4 * n), ncol = n), 1, stats::median)
    q <- stats::quantile(medians, c(0.95, 0.99), names = FALSE)
    expect_true(
      all(abs(.median_normal_sf(q, n) - c(0.05, 0.01)) < 0.004),
      info = paste("n =", n)
    )
  }
})

test_that(".peptide_z maps t-statistics to standard normal scores", {
  t <- c(-4, -1, 0, 2)
  expect_equal(.peptide_z(t, 6), stats::qnorm(stats::pt(t, 6)))
  expect_equal(.peptide_z(1.5, Inf), 1.5)
  expect_true(is.finite(.peptide_z(80, 6)))
})

test_that(".ropeca_v2_pvalue is uniform under H0 for independent and correlated peptides", {
  for (rho in c(0, 0.3)) {
    for (n in c(2, 3, 5, 10, 20)) {
      median_z <- simulate_median_z(4000, n, rho, seed = n)
      p <- .ropeca_v2_pvalue(median_z, rep(n, length(median_z)), rho)
      expect_true(abs(mean(p < 0.05) - 0.05) < 0.015, info = paste("rho =", rho, "n =", n))
    }
  }
})

test_that(".ropeca_v2_pvalue reduces to the peptide test and to the independent null", {
  z <- c(-3, -1, 0, 0.5, 2.5)
  expect_equal(.ropeca_v2_pvalue(z, rep(1, 5), 0.3), 2 * stats::pnorm(-abs(z)), tolerance = 1e-6)
  m <- c(0.2, 1, 2)
  expect_equal(.ropeca_v2_pvalue(m, rep(4, 3), 0), pmin(1, 2 * .median_normal_sf(m, 4)), tolerance = 1e-8)
})

test_that(".ropeca_v2_rho estimates the within-protein peptide correlation", {
  muffle_warnings_matching(one_group_warning, expr = {
    strategy <- strategy_lm("abundance ~ group_")
    subject_id <- c("protein_Id", "peptide_Id")

    correlated <- simulate_correlated_peptides(300, 5, rho = 0.3, missing = 0.3)
    model <- build_model(correlated, strategy, subject_id = subject_id)
    expect_true(abs(.ropeca_v2_rho(model$model_df, subject_id, "sampleName") - 0.3) < 0.05)

    independent <- simulate_correlated_peptides(300, 5, rho = 0, seed = 2)
    model <- build_model(independent, strategy, subject_id = subject_id)
    expect_true(.ropeca_v2_rho(model$model_df, subject_id, "sampleName") < 0.03)
  })
})

test_that("ContrastsROPECAV2 checks its arguments", {
  peptide_contrast <- list(subject_id = c("protein_Id", "peptide_Id"))
  expect_error(ContrastsROPECAV2$new(peptide_contrast, rho = 1), class = "prolfqua_error_bad_argument")
  expect_error(ContrastsROPECAV2$new(peptide_contrast, rho = NA_real_), class = "prolfqua_error_bad_argument")
  expect_error(
    ContrastsROPECAV2$new(list(subject_id = "protein_Id"), rho = 0.2),
    class = "prolfqua_error_bad_argument"
  )
})
