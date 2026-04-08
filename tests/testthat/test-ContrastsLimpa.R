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
  expect_equal(length(lfq_agg$config$hierarchy_keys()), 1)
  expect_equal(lfq_agg$config$hierarchy_keys(), "protein_Id")

  # Should have no NAs in intensity
  expect_equal(sum(is.na(lfq_agg$data[[lfq_agg$config$get_response()]])), 0)

  # Response column should be "limpa"
  expect_equal(lfq_agg$config$get_response(), "limpa")

  # Should have SE column
  se_col <- lfq_agg$config$opt_se
  expect_true(nchar(se_col) > 0)
  expect_true(se_col %in% colnames(lfq_agg$data))
  expect_equal(sum(is.na(lfq_agg$data[[se_col]])), 0)

  # Should have nr_children (observation counts)
  nr_col <- lfq_agg$config$nr_children
  expect_true(nchar(nr_col) > 0)
  expect_true(nr_col %in% colnames(lfq_agg$data))

  # DPC should be estimated
  expect_true(!is.null(agg$dpc_result))

  # Should be pivotable to wide
  wide <- lfq_agg$to_wide(as.matrix = TRUE)
  expect_true(is.matrix(wide$data))
  expect_equal(sum(is.na(wide$data)), 0)

  # SE should also be pivotable
  se_wide <- lfq_agg$to_wide(as.matrix = TRUE, value = se_col)
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
    length(lfq_imp$config$hierarchy_keys()),
    length(lfqdata$config$hierarchy_keys())
  )

  # Should have no NAs in intensity
  expect_equal(sum(is.na(lfq_imp$data[[lfq_imp$config$get_response()]])), 0)

  # Should have SE column
  expect_true(nchar(lfq_imp$config$opt_se) > 0)
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
