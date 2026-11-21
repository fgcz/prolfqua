# Convergence controls: settable per strategy and aggregation method, documented defaults, classed errors on
# unknown names.

#' @srrstats {RE3.3} Tests set iteration limits and tolerances for rlm, lmer, binomial, Firth, logistf, rfit and
#'   the medpolish and rlm aggregation, and check that unknown control names are rejected.
#' @noRd
NULL

conv_protein <- function() {
  sim <- sim_lfq_data_peptide_config(n_proteins = 10, with_missing = FALSE, seed = 42)
  peptide <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  list(peptide = peptide, protein = suppressWarnings(peptide$get_Aggregator()$aggregate()))
}

test_that("strategies keep their documented defaults and merge user controls", {
  expect_equal(strategy_rlm("y ~ g")$control, list(maxit = 20, acc = 1e-4))
  expect_equal(strategy_rlm("y ~ g", control = list(maxit = 5))$control, list(maxit = 5, acc = 1e-4))
  expect_equal(strategy_binomial("~ g")$control, list(epsilon = 1e-8, maxit = 25))
  expect_equal(strategy_binomial_firth("~ g")$control, list(maxit = 100, tol = 1e-8))
  expect_equal(strategy_rfit("y ~ g")$control, list(maxit = 200, reltol = .Machine$double.eps^(3 / 4)))
  expect_equal(strategy_lmer("y ~ g + (1 | p)")$control, list())
  expect_equal(strategy_lm("y ~ g")$control, list())
})

test_that("unknown or unnamed controls raise a classed error", {
  expect_error(strategy_rlm("y ~ g", control = list(maxiter = 5)), "maxiter", class = "prolfqua_error_bad_argument")
  expect_error(strategy_rlm("y ~ g", control = list(5)), class = "prolfqua_error_bad_argument")
  expect_error(strategy_binomial_firth("~ g", control = list(epsilon = 1)), class = "prolfqua_error_bad_argument")
  expect_error(StrategyBase$new("y ~ g", control = list(maxit = 1)), "nothing", class = "prolfqua_error_bad_argument")
  skip_if_not_installed("lme4")
  expect_error(strategy_lmer("y ~ g + (1 | p)", control = list(nonsense = 1)), class = "prolfqua_error_bad_argument")
})

test_that("an rlm iteration limit reaches MASS::rlm through build_contrast_analysis", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- conv_protein()$protein
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  facade <- suppressWarnings(build_contrast_analysis(protein, "~ group_", contrasts, "rlm", control = list(maxit = 1)))
  fits <- facade$model$model_df$linear_model[facade$model$model_df$has_model_fit]
  expect_true(all(vapply(fits, function(f) f$control$maxit %||% length(f$conv), numeric(1)) <= 1))
  failures <- facade$model$provenance()
  expect_true(any(grepl("'rlm' failed to converge in 1 steps", failures$reason)))
})

test_that("binomial, Firth and rfit honour their controls", {
  dat <- data.frame(
    group_ = factor(rep(c("A", "B"), each = 4)),
    detected = c(1, 2, 1, 3, 4, 5, 3, 5),
    undetected = c(4, 3, 4, 2, 1, 0, 2, 0)
  )
  glm_fit <- suppressWarnings(strategy_binomial("~ group_", control = list(maxit = 1))$model_fun(dat))
  expect_equal(glm_fit$control$maxit, 1)
  expect_match(strategy_binomial_firth("~ group_", control = list(maxit = 1))$model_fun(dat), "did not converge")
  firth <- strategy_binomial_firth("~ group_")$model_fun(dat)
  expect_gt(firth$firth_iterations, 1)
  skip_if_not_installed("Rfit")
  rfit_data <- data.frame(y = c(1.2, 0.8, 1.1, 0.9, 2.1, 1.8, 2.3, 2.0), g = factor(rep(c("a", "b"), each = 4)))
  expect_warning(strategy_rfit("y ~ g", control = list(maxit = 1))$model_fun(rfit_data), "Convergence status")
})

test_that("aggregation iteration limits reach medpolish and rlm", {
  withr::local_options(prolfqua.verbose = FALSE)
  peptide <- conv_protein()$peptide
  expect_error(peptide$get_Aggregator("rlm", control = list(maxiter = 1)), class = "prolfqua_error_bad_argument")
  for (method in c("medpolish", "rlm")) {
    control <- if (method == "medpolish") list(maxiter = 1, eps = 0) else list(maxit = 1)
    messages <- character(0)
    withCallingHandlers(
      peptide$get_Aggregator(method, control = control)$aggregate(),
      prolfqua_fit_failed = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_length(messages, 1L)
    expect_match(messages, "converge in 1 ")
  }
})

test_that("default controls leave aggregation and rlm results unchanged", {
  withr::local_options(prolfqua.verbose = FALSE)
  peptide <- conv_protein()$peptide
  default <- suppressWarnings(peptide$get_Aggregator("rlm")$aggregate())
  explicit <- suppressWarnings(peptide$get_Aggregator("rlm", control = list(maxit = 20, acc = 1e-4))$aggregate())
  expect_equal(default$data_long(), explicit$data_long())
})

#' @srrstats {RE4.7} Tests read per-subject convergence status and iteration counts from the model table and the
#'   not_converged count from the fit status.
#' @noRd
NULL

test_that("model tables record convergence and iterations per subject", {
  withr::local_options(prolfqua.verbose = FALSE)
  fixture <- conv_protein()
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  lm <- build_contrast_analysis(fixture$protein, "~ group_", contrasts, "lm")
  expect_true(all(lm$model$model_df$converged))
  expect_true(all(is.na(lm$model$model_df$iterations)))
  rlm <- suppressWarnings(build_contrast_analysis(
    fixture$protein,
    "~ group_",
    contrasts,
    "rlm",
    control = list(maxit = 1)
  ))
  model_df <- rlm$model$model_df
  expect_type(model_df$converged, "logical")
  expect_type(model_df$iterations, "integer")
  expect_true(all(model_df$iterations[model_df$has_model_fit] <= 1L))
  status <- summary(rlm)$fit_status
  expect_equal(status$not_converged, sum(!model_df$converged, na.rm = TRUE))
  expect_gt(status$not_converged, 0L)
  default <- suppressWarnings(build_contrast_analysis(fixture$protein, "~ group_", contrasts, "rlm"))
  expect_true(all(default$model$model_df$iterations >= 1L, na.rm = TRUE))
})

test_that("binomial and Firth detection models record glm and Firth iterations", {
  withr::local_options(prolfqua.verbose = FALSE)
  sim <- sim_lfq_data_peptide_config(n_proteins = 10, weight_missing = 0.4, seed = 42)
  peptide <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  for (method in c("binomial_nested", "binomial_firth_nested")) {
    facade <- suppressWarnings(build_contrast_analysis(peptide, "~ group_", contrasts, method))
    model_df <- facade$model$model_df
    fitted <- model_df$has_model_fit
    expect_true(all(model_df$converged[fitted]), label = method)
    expect_true(all(model_df$iterations[fitted] >= 1L), label = method)
  }
})

test_that("lmer fits record the optimizer status", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  withr::local_options(prolfqua.verbose = FALSE)
  peptide <- conv_protein()$peptide
  facade <- suppressMessages(suppressWarnings(build_contrast_analysis(
    peptide,
    "~ group_ + (1 | peptide_Id)",
    c(AvsC = "group_A - group_Ctrl"),
    "lmer_nested"
  )))
  model_df <- facade$model$model_df
  fitted <- model_df$has_model_fit
  expect_type(model_df$converged, "logical")
  expect_true(all(model_df$iterations[fitted] > 0L))
})
