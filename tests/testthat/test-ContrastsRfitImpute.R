## Tests for the rfit_impute backend
## (new_imputed_model on an rfit fit, build_model_impute with strategy_rfit,
##  ContrastsRfitImputeFacade, build_contrast_analysis dispatch, registry)

make_protein_lfqdata_rfit_impute <- function(Nprot = 50, weight_missing = 0.5, seed = 42) {
  istar <- prolfqua::sim_lfq_data_protein_config(
    Nprot = Nprot,
    weight_missing = weight_missing,
    seed = seed
  )
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  lfqdata
}

RFITIMP_CONTRASTS <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
RFITIMP_MODELSTR <- "~ group_"


test_that("new_imputed_model wraps an rfit fit, overriding vcov/sigma/df.residual", {
  skip_if_not_installed("Rfit")
  set.seed(1)
  d <- data.frame(
    y = rnorm(12),
    group_ = factor(rep(c("A", "Ctrl"), each = 6))
  )
  fit <- prolfqua::strategy_rfit("y ~ group_")$model_fun(d)
  expect_true(inherits(fit, "rfit_prolfqua"))

  borrowed_vcov <- stats::vcov(fit) * 2
  wrapped <- prolfqua:::new_imputed_model(
    fit,
    borrowed_vcov = borrowed_vcov,
    borrowed_sigma = 0.99,
    borrowed_df = 5,
    n_observed = 8
  )

  # Backend-neutral class is prepended; rfit behaviour is preserved underneath.
  expect_s3_class(wrapped, "imputed_model")
  expect_s3_class(wrapped, "rfit_prolfqua")

  # vcov / sigma / df.residual dispatch to the borrowed values.
  expect_equal(stats::vcov(wrapped), borrowed_vcov)
  expect_equal(stats::sigma(wrapped), 0.99)
  expect_equal(stats::df.residual(wrapped), 5)

  # coef / terms fall through to the wrapped rfit fit.
  expect_equal(stats::coef(wrapped), stats::coef(fit))
  expect_equal(stats::terms(wrapped), stats::terms(fit))
})


test_that("build_model_impute rescues failed/singular rfit fits (vcov, fail-hard)", {
  skip_if_not_installed("Rfit")
  lfqdata <- make_protein_lfqdata_rfit_impute()
  strat <- prolfqua::strategy_rfit(paste(lfqdata$response(), RFITIMP_MODELSTR))

  mod_plain <- prolfqua::build_model(lfqdata, strat)
  mod_impute <- prolfqua::build_model_impute(
    lfqdata,
    strat,
    borrow_method = "vcov",
    on_misalign = "fail"
  )

  # The rescue adds an `imputed` flag and recovers at least some proteins.
  expect_true("imputed" %in% colnames(mod_impute$model_df))
  n_imputed <- sum(mod_impute$model_df$imputed, na.rm = TRUE)
  expect_gt(n_imputed, 0)
  expect_false("imputed" %in% colnames(mod_plain$model_df))

  n_good_plain <- sum(
    mod_plain$model_df$has_model_fit &
      !is.na(mod_plain$model_df$isSingular) &
      !mod_plain$model_df$isSingular
  )
  n_good_impute <- sum(
    mod_impute$model_df$has_model_fit &
      !is.na(mod_impute$model_df$isSingular) &
      !mod_impute$model_df$isSingular
  )
  expect_gte(n_good_impute, n_good_plain)
})


test_that("ContrastsRfitImputeFacade returns the standard schema with facade tag", {
  skip_if_not_installed("Rfit")
  lfqdata <- make_protein_lfqdata_rfit_impute()
  fa <- prolfqua::ContrastsRfitImputeFacade$new(lfqdata, RFITIMP_MODELSTR, RFITIMP_CONTRASTS)

  expect_true(inherits(fa, "ContrastsRfitImputeFacade"))
  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  for (col in c(
    "modelName",
    "estimate_type",
    "contrast",
    "diff",
    "FDR",
    "p.value",
    "statistic",
    "std.error",
    "df",
    "conf.low",
    "conf.high",
    "sigma"
  )) {
    expect_true(col %in% colnames(res), info = paste("missing column:", col))
  }
  expect_false("facade" %in% colnames(res))
  expect_true(!any(is.na(res$diff)))
  expect_true(all(res$modelName == "rfit_impute"))

  # Interface methods all work.
  expect_true(is.data.frame(fa$get_missing()))
  expect_true(nrow(fa$to_wide()) > 0)
  expect_s3_class(fa$get_Plotter(), "ContrastsPlotter")
})


test_that("rfit_impute flags rescued rows via estimate_type", {
  skip_if_not_installed("Rfit")
  lfqdata <- make_protein_lfqdata_rfit_impute()
  fa <- prolfqua::ContrastsRfitImputeFacade$new(lfqdata, RFITIMP_MODELSTR, RFITIMP_CONTRASTS)
  res <- fa$get_contrasts()
  # modelName is the facade key; rescued rows are tagged in estimate_type.
  expect_setequal(unique(res$modelName), "rfit_impute")
  expect_true("observed" %in% res$estimate_type)
  expect_true("lod_imputed" %in% res$estimate_type)
})


test_that("rfit_impute has fewer or equal missing than plain rfit", {
  skip_if_not_installed("Rfit")
  lfqdata <- make_protein_lfqdata_rfit_impute()
  fa_impute <- prolfqua::ContrastsRfitImputeFacade$new(lfqdata, RFITIMP_MODELSTR, RFITIMP_CONTRASTS)
  fa_rfit <- prolfqua::ContrastsRfitFacade$new(lfqdata, RFITIMP_MODELSTR, RFITIMP_CONTRASTS)

  expect_gte(nrow(fa_impute$get_contrasts()), nrow(fa_rfit$get_contrasts()))
  expect_lte(nrow(fa_impute$get_missing()), nrow(fa_rfit$get_missing()))
})


test_that("rfit_impute is stable under LOD-created ties", {
  skip_if_not_installed("Rfit")
  # High missingness creates whole missing groups; LOD imputation then yields
  # many identical (tied) values whose rank order Rfit resolves with
  # ties.method = "first". Deterministic sample-row ordering must make the
  # facade reproducible across runs.
  lfqdata <- make_protein_lfqdata_rfit_impute(Nprot = 40, weight_missing = 1, seed = 7)
  res1 <- prolfqua::ContrastsRfitImputeFacade$new(
    lfqdata,
    RFITIMP_MODELSTR,
    RFITIMP_CONTRASTS
  )$get_contrasts()
  res2 <- prolfqua::ContrastsRfitImputeFacade$new(
    lfqdata,
    RFITIMP_MODELSTR,
    RFITIMP_CONTRASTS
  )$get_contrasts()
  expect_equal(res1, res2)
})


test_that("rfit null-model fallback names survive the imputed wrapper", {
  skip_if_not_installed("Rfit")
  # Two groups whose medians do not separate -> Rfit returns the unnamed
  # intercept-only fallback; the strategy restores names from the design
  # matrix (rfit.R). That fix must hold for refits that get wrapped too.
  d <- data.frame(
    y = c(13, 15, 17, 14, 15, 16),
    group_ = factor(rep(c("A", "Ctrl"), each = 3))
  )
  # Factor levels c("A", "Ctrl") -> reference "A" -> coefficient "group_Ctrl".
  fit <- prolfqua::strategy_rfit("y ~ group_")$model_fun(d)
  expect_equal(names(stats::coef(fit)), c("(Intercept)", "group_Ctrl"))

  wrapped <- prolfqua:::new_imputed_model(
    fit,
    borrowed_vcov = stats::vcov(fit),
    borrowed_sigma = stats::sigma(fit),
    borrowed_df = 4,
    n_observed = 4
  )
  expect_equal(names(stats::coef(wrapped)), c("(Intercept)", "group_Ctrl"))
  lf <- prolfqua::linfct_from_model(wrapped, as_list = FALSE)
  expect_false(is.null(colnames(lf)))
  expect_setequal(colnames(lf), c("(Intercept)", "group_Ctrl"))
})


test_that("build_contrast_analysis dispatches method = 'rfit_impute'", {
  skip_if_not_installed("Rfit")
  lfqdata <- make_protein_lfqdata_rfit_impute()
  fa <- prolfqua::build_contrast_analysis(
    lfqdata,
    RFITIMP_MODELSTR,
    RFITIMP_CONTRASTS,
    method = "rfit_impute"
  )
  expect_true(inherits(fa, "ContrastsRfitImputeFacade"))
  expect_true(all(fa$get_contrasts()$modelName == "rfit_impute"))
})


test_that("rfit_impute facade is registered as needs = 'same'", {
  entry <- prolfqua::lookup_facade("rfit_impute")
  expect_equal(entry$class, "ContrastsRfitImputeFacade")
  expect_equal(entry$needs, "same")
  expect_true("rfit_impute" %in% names(prolfqua::list_facades()))
  expect_equal(
    prolfqua::FACADE_REGISTRY$rfit_impute$class,
    "ContrastsRfitImputeFacade"
  )
})
