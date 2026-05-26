test_that(".call_limpa_dpc supports limpa versions with and without b1.upper", {
  expr_matrix <- matrix(seq_len(6), nrow = 2)

  dpc_without_b1_upper <- function(y) {
    list(y = y)
  }
  expect_equal(
    prolfqua:::.call_limpa_dpc(expr_matrix, dpc_fun = dpc_without_b1_upper)$y,
    expr_matrix
  )

  dpc_with_b1_upper <- function(y, b1.upper = NULL) {
    list(y = y, b1.upper = b1.upper)
  }
  expect_equal(
    prolfqua:::.call_limpa_dpc(expr_matrix, b1.upper = 1, dpc_fun = dpc_with_b1_upper)$b1.upper,
    1
  )
})


test_that("strategy_limpa returns StrategyLimpa R6 object", {
  skip_if_not_installed("limpa")
  strat <- prolfqua::strategy_limpa("abundance ~ group_")
  expect_true(inherits(strat, "StrategyLimpa"))
  expect_true(inherits(strat$formula, "formula"))
  expect_equal(strat$model_name, "limpa")
  expect_false(strat$trend)
  expect_false(strat$robust)
})


test_that("AggregateLimpa produces protein-level LFQData with SE and n_obs", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  agg <- prolfqua::AggregateLimpa$new(lfqdata, "protein")
  lfq_agg <- agg$aggregate()

  # Should be protein-level (reduced hierarchy)
  expect_equal(length(lfq_agg$get_config()$hierarchy_keys()), 1)
  expect_equal(lfq_agg$get_config()$hierarchy_keys(), "protein_Id")

  # Should have no NAs in intensity
  expect_equal(sum(is.na(lfq_agg$data_long()[[lfq_agg$get_config()$get_response()]])), 0)

  # Response column should be "limpa"
  expect_equal(lfq_agg$get_config()$get_response(), "limpa")

  # Should have SE column
  se_col <- lfq_agg$get_config()$opt_se
  expect_true(nchar(se_col) > 0)
  expect_true(se_col %in% colnames(lfq_agg$data_long()))
  expect_equal(sum(is.na(lfq_agg$data_long()[[se_col]])), 0)

  # Should have nr_children (observation counts)
  nr_col <- lfq_agg$get_config()$nr_children
  expect_true(nchar(nr_col) > 0)
  expect_true(nr_col %in% colnames(lfq_agg$data_long()))

  # DPC should be estimated
  expect_true(!is.null(agg$dpc_result))

  # Should be pivotable to wide
  wide <- lfq_agg$data_wide(as.matrix = TRUE)
  expect_true(is.matrix(wide$data))
  expect_equal(sum(is.na(wide$data)), 0)

  # SE should also be pivotable
  se_wide <- lfq_agg$data_wide(as.matrix = TRUE, value = se_col)
  expect_true(is.matrix(se_wide$data))
  expect_equal(dim(se_wide$data), dim(wide$data))
})


test_that("AggregateLimpa impute_only mode preserves hierarchy", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  agg <- prolfqua::AggregateLimpa$new(lfqdata, impute_only = TRUE)
  lfq_imp <- agg$aggregate()

  # Hierarchy should be the same as input (not reduced)
  expect_equal(
    length(lfq_imp$get_config()$hierarchy_keys()),
    length(lfqdata$get_config()$hierarchy_keys())
  )

  # Should have no NAs in intensity
  expect_equal(sum(is.na(lfq_imp$data_long()[[lfq_imp$get_config()$get_response()]])), 0)

  # Should have SE column
  expect_true(nchar(lfq_imp$get_config()$opt_se) > 0)
})


test_that("build_model_limpa returns ModelLimma", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  agg <- prolfqua::AggregateLimpa$new(lfqdata, "protein")
  lfq_agg <- agg$aggregate()

  strat <- prolfqua::strategy_limpa("limpa ~ group_")
  mod <- prolfqua::build_model_limpa(lfq_agg, strat)

  expect_true(inherits(mod, "ModelLimma"))
  expect_true(inherits(mod, "ModelInterface"))
  expect_true(!is.null(mod$fit))
  expect_true(!is.null(mod$dummy_model))
  expect_true(!is.null(mod$rowdata))
  expect_equal(mod$model_name, "limpa")

  # Coefficients should work
  coeffs <- mod$get_coefficients()
  expect_true(is.data.frame(coeffs))
  expect_true(nrow(coeffs) > 0)
})


test_that("ContrastsLimma works with limpa ModelLimma", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  agg <- prolfqua::AggregateLimpa$new(lfqdata, "protein")
  lfq_agg <- agg$aggregate()

  strat <- prolfqua::strategy_limpa("limpa ~ group_")
  mod <- prolfqua::build_model_limpa(lfq_agg, strat)

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  con <- prolfqua::ContrastsLimma$new(mod, contrasts)

  res <- con$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true("diff" %in% colnames(res))
  expect_true("p.value" %in% colnames(res))
  expect_true("statistic" %in% colnames(res))
  expect_true(nrow(res) > 0)
})


test_that("ContrastsLimpaFacade end-to-end", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  agg <- prolfqua::AggregateLimpa$new(lfqdata, "protein")
  lfq_agg <- agg$aggregate()

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::ContrastsLimpaFacade$new(lfq_agg, "~ group_", contrasts)

  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true("diff" %in% colnames(res))
  expect_true("p.value" %in% colnames(res))
  expect_true("FDR" %in% colnames(res))
  expect_true("modelName" %in% colnames(res))
  expect_true(all(res$modelName == "limpa"))
  expect_true(nrow(res) > 0)

  # get_missing should work
  miss <- fa$get_missing()
  expect_true(is.data.frame(miss))

  # to_wide should work
  wide <- fa$to_wide()
  expect_true(is.data.frame(wide))
})


test_that("build_contrast_analysis with method='limpa' works", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  agg <- prolfqua::AggregateLimpa$new(lfqdata, "protein")
  lfq_agg <- agg$aggregate()

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::build_contrast_analysis(lfq_agg, "~ group_", contrasts, method = "limpa")

  expect_true(inherits(fa, "ContrastsLimpaFacade"))
  res <- fa$get_contrasts()
  expect_true(nrow(res) > 0)
  expect_true("diff" %in% colnames(res))
})


test_that("ContrastsLimpaNestedFacade end-to-end (precursor input)", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::ContrastsLimpaNestedFacade$new(lfqdata, "~ group_", contrasts)

  expect_true(inherits(fa, "ContrastsLimpaNestedFacade"))
  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  expect_true(all(res$facade == "limpa_nested"))
  expect_true("diff" %in% colnames(res))
  expect_true("p.value" %in% colnames(res))
  expect_true("FDR" %in% colnames(res))
})


test_that("build_contrast_analysis with method='limpa_nested' works", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limpa_nested")

  expect_true(inherits(fa, "ContrastsLimpaNestedFacade"))
  res <- fa$get_contrasts()
  expect_true(nrow(res) > 0)
})


test_that("ContrastsLimpaFacade rejects LFQData without opt_se", {
  skip_if_not_installed("limpa")
  # Use protein-level data without SE column
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  expect_error(
    prolfqua::ContrastsLimpaFacade$new(lfqdata, "~ group_", contrasts),
    "opt_se"
  )
})
