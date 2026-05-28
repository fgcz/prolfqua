## Tests for the Rfit rank-based regression backend
## (StrategyRfit, S3 glue, ContrastsRfitFacade, build_contrast_analysis dispatch)

make_protein_lfqdata_rfit <- function(Nprot = 30) {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot, seed = 42)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  lfqdata
}

RFIT_CONTRASTS <- c("A_vs_Ctrl" = "group_A - group_Ctrl", "B_vs_Ctrl" = "group_B - group_Ctrl")
RFIT_MODELSTR <- "~ group_"


test_that("strategy_rfit constructs a StrategyRfit", {
  strat <- prolfqua::strategy_rfit("Intensity ~ group_")
  expect_true(inherits(strat, "StrategyRfit"))
  expect_false(strat$is_mixed)
  expect_equal(strat$model_fun(get_formula = TRUE), stats::as.formula("Intensity ~ group_"))
})


test_that("rfit S3 glue: vcov is named, df.residual and sigma resolve", {
  skip_if_not_installed("Rfit")
  set.seed(1)
  d <- data.frame(
    y = rnorm(24),
    Treatment = factor(rep(c("A", "B"), each = 12)),
    Background = factor(rep(c("X", "Z"), times = 12))
  )
  fo <- y ~ Treatment * Background
  strat <- prolfqua::strategy_rfit("y ~ Treatment * Background")
  fit <- strat$model_fun(d)

  expect_true(inherits(fit, "rfit_prolfqua"))
  v <- stats::vcov(fit)
  expect_equal(rownames(v), names(stats::coef(fit)))
  expect_equal(colnames(v), names(stats::coef(fit)))
  expect_equal(stats::df.residual(fit), nrow(d) - length(stats::coef(fit)))
  expect_equal(stats::sigma(fit), fit$tauhat)
})


test_that("linfct_from_model on an rfit fit matches the lm scaffold", {
  skip_if_not_installed("Rfit")
  set.seed(2)
  d <- data.frame(
    y = rnorm(24),
    Treatment = factor(rep(c("A", "B"), each = 12)),
    Background = factor(rep(c("X", "Z"), times = 12))
  )
  fo <- y ~ Treatment * Background
  fit <- prolfqua::strategy_rfit("y ~ Treatment * Background")$model_fun(d)
  m <- lm(fo, data = d)

  lf_rfit <- prolfqua::linfct_from_model(fit, as_list = FALSE)
  lf_lm <- prolfqua::linfct_from_model(m, as_list = FALSE)
  expect_equal(lf_rfit, lf_lm)
})


test_that("ContrastsRfitFacade initialises and returns the standard schema", {
  skip_if_not_installed("Rfit")
  lfqdata <- make_protein_lfqdata_rfit()
  fa <- prolfqua::ContrastsRfitFacade$new(lfqdata, RFIT_MODELSTR, RFIT_CONTRASTS)

  expect_true(inherits(fa, "ContrastsRfitFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))

  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  for (col in c("facade", "modelName", "contrast", "diff", "FDR", "p.value", "statistic", "std.error")) {
    expect_true(col %in% colnames(res), info = paste("missing column:", col))
  }
  expect_true(!any(is.na(res$diff)))
  expect_true(!any(is.na(res$statistic)))
  expect_true(all(res$facade == "rfit"))

  expect_true(is.data.frame(fa$get_missing()))
  expect_true(nrow(fa$to_wide()) > 0)
  expect_true(!is.null(fa$get_Plotter()))
})


test_that("rfit fold-changes agree in sign and correlate with lm", {
  skip_if_not_installed("Rfit")
  lfqdata <- make_protein_lfqdata_rfit()
  res_rfit <- prolfqua::ContrastsRfitFacade$new(lfqdata, RFIT_MODELSTR, RFIT_CONTRASTS)$get_contrasts()
  res_lm <- prolfqua::ContrastsLMFacade$new(lfqdata, RFIT_MODELSTR, RFIT_CONTRASTS)$get_contrasts()

  mrg <- merge(
    res_rfit[, c("protein_Id", "contrast", "diff")],
    res_lm[, c("protein_Id", "contrast", "diff")],
    by = c("protein_Id", "contrast"),
    suffixes = c(".rfit", ".lm")
  )
  expect_true(nrow(mrg) > 0)
  expect_gt(mean(sign(mrg$diff.rfit) == sign(mrg$diff.lm)), 0.9)
  expect_gt(cor(mrg$diff.rfit, mrg$diff.lm), 0.9)
})


test_that("build_contrast_analysis dispatches method = 'rfit'", {
  skip_if_not_installed("Rfit")
  lfqdata <- make_protein_lfqdata_rfit()
  fa <- prolfqua::build_contrast_analysis(lfqdata, RFIT_MODELSTR, RFIT_CONTRASTS, method = "rfit")
  expect_true(inherits(fa, "ContrastsRfitFacade"))
  expect_true(all(fa$get_contrasts()$facade == "rfit"))
})


test_that("rfit facade is registered as needs = 'same'", {
  entry <- prolfqua::lookup_facade("rfit")
  expect_equal(entry$class, "ContrastsRfitFacade")
  expect_equal(entry$needs, "same")
})


test_that("rfit handles a two-factor design with interactions", {
  skip_if_not_installed("Rfit")
  istar <- prolfqua::sim_lfq_data_2factor_config(Nprot = 12, with_missing = FALSE)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  contrasts <- c("TB_vs_A" = "TreatmentB - TreatmentA", "BgZ_vs_X" = "BackgroundZ - BackgroundX")
  fa <- prolfqua::ContrastsRfitFacade$new(lfqdata, "~ Treatment * Background", contrasts)
  res <- fa$get_contrasts()
  expect_setequal(unique(res$contrast), names(contrasts))
  expect_true(!any(is.na(res$statistic)))
})
