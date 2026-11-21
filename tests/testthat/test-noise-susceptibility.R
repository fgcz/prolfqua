# Noise susceptibility: trivial input noise and different random seeds must not change results.

#' @srrstats {G5.9} Noise susceptibility is tested for one method per fitting engine (lm, rlm, limma, lmer_nested,
#'   firth) in the standard suite and for every registered modelling method in the extended suite.
#' @srrstats {G5.9a} Multiplying the response by 1 + .Machine$double.eps * N(0, 1) changes contrast estimates,
#'   standard errors and p-values by less than 1e-6 (relative); 1e-4 for limpa, whose iterative estimation of
#'   protein quantities and weights moves at about the 1e-5 level.
#' @srrstats {G5.9b} Running the same analysis after different set.seed() calls gives identical results; the
#'   fitting algorithms contain no random component.
#' @noRd
NULL

NOISE_WARNINGS <- paste(
  rlm_convergence_warning,
  loess_ties_warning,
  single_peptide_warning,
  one_group_warning,
  partial_na_warning,
  "Zero sample variances detected",
  sep = "|"
)
# limpa estimates protein quantities and weights iteratively; trivial input noise moves its results at ~1e-5.
NOISE_TOLERANCE <- c(limpa = 1e-4, limpa_nested = 1e-4)
# Standard suite: one method per fitting engine; the extended suite runs all registered methods.
NOISE_STANDARD_METHODS <- c("lm", "rlm", "limma", "lmer_nested", "firth")
NOISE_PACKAGES <- list(
  rfit = "Rfit",
  rfit_impute = "Rfit",
  lmer_nested = c("lme4", "lmerTest"),
  firth = "logistf",
  limpa = "limpa",
  limpa_nested = "limpa"
)

noise_case_input <- function(method, perturb = FALSE) {
  sim <- sim_lfq_data_peptide_config(n_proteins = 20, seed = 42)
  data <- sim$data
  if (perturb) {
    response <- sim$config$get_response()
    noise <- withr::with_seed(7, stats::rnorm(nrow(data)))
    data[[response]] <- data[[response]] * (1 + .Machine$double.eps * noise)
  }
  peptide <- LFQData$new(data, sim$config)$get_Transformer()$log2()$lfq
  if (method == "limpa") {
    return(AggregateLimpa$new(peptide, "protein")$aggregate())
  }
  if (lookup_facade(method)$needs == "nested") {
    return(peptide)
  }
  peptide$get_Aggregator()$aggregate()
}

noise_case_contrasts <- function(method, perturb = FALSE) {
  modelstr <- if (method == "lmer_nested") "~ group_ + (1 | peptide_Id) + (1 | sampleName)" else "~ group_"
  muffle_warnings_matching(
    {
      res <- build_contrast_analysis(
        noise_case_input(method, perturb),
        modelstr,
        c(A_vs_Ctrl = "group_A - group_Ctrl"),
        method
      )$get_contrasts()
      res[order(res$protein_Id), , drop = FALSE]
    },
    NOISE_WARNINGS
  )
}

numeric_result_columns <- function(res) {
  intersect(c("diff", "std.error", "statistic", "p.value"), names(res))
}

for (method in names(list_facades())) {
  test_that(paste("trivial noise and seeds leave", method, "results unchanged"), {
    if (!method %in% NOISE_STANDARD_METHODS) {
      skip_if_not_extended()
    }
    for (package in NOISE_PACKAGES[[method]]) {
      skip_if_not_installed(package)
    }
    original <- noise_case_contrasts(method)
    perturbed <- noise_case_contrasts(method, perturb = TRUE)
    columns <- numeric_result_columns(original)
    expect_gt(length(columns), 0L)
    expect_identical(perturbed$protein_Id, original$protein_Id)
    tolerance <- if (method %in% names(NOISE_TOLERANCE)) NOISE_TOLERANCE[[method]] else 1e-6
    for (column in columns) {
      expect_equal(perturbed[[column]], original[[column]], tolerance = tolerance, label = paste(method, column))
    }
    first <- withr::with_seed(1, noise_case_contrasts(method))
    second <- withr::with_seed(2, noise_case_contrasts(method))
    expect_identical(first[columns], second[columns])
  })
}
