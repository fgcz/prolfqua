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
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 20)
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

  # DPC should be estimated: intercept and slope of the detection probability curve, from all samples
  expect_named(agg$dpc_result$dpc, c("beta0", "beta1"))
  expect_equal(agg$dpc_result$nsamples, dplyr::n_distinct(lfqdata$data_long()[[lfqdata$sample_name()]]))

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
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 20)
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
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  agg <- prolfqua::AggregateLimpa$new(lfqdata, "protein")
  lfq_agg <- agg$aggregate()

  strat <- prolfqua::strategy_limpa("limpa ~ group_")
  mod <- prolfqua::build_model_limpa(lfq_agg, strat)

  expect_true(inherits(mod, "ModelLimma"))
  expect_true(inherits(mod, "ModelInterface"))
  # one fit row per protein, one coefficient per level of group_
  expect_s4_class(mod$fit, "MArrayLM")
  expect_equal(dim(mod$fit$coefficients), c(20L, 3L))
  expect_equal(names(stats::coef(mod$dummy_model)), colnames(mod$fit$coefficients))
  expect_equal(mod$rowdata$protein_Id, rownames(mod$fit$coefficients))
  expect_equal(mod$model_name, "limpa")

  # Coefficients should work: 20 proteins x 3 coefficients
  coeffs <- mod$get_coefficients()
  expect_true(is.data.frame(coeffs))
  expect_equal(nrow(coeffs), 20 * 3)
})


test_that("ContrastsLimma works with limpa ModelLimma", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 30)
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
  # limpa imputes every protein, so each one has a contrast
  expect_equal(nrow(res), dplyr::n_distinct(lfq_agg$data_long()$protein_Id))
})


test_that("ContrastsLimpaFacade end-to-end", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 30)
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
  n_proteins <- dplyr::n_distinct(lfq_agg$data_long()$protein_Id)
  expect_equal(nrow(res), n_proteins)

  # get_missing should work
  miss <- fa$get_missing()
  expect_true(is.data.frame(miss))

  # to_wide should work
  wide <- fa$to_wide()
  expect_true(is.data.frame(wide))
  expect_equal(nrow(wide), n_proteins)
})


test_that("build_contrast_analysis with method='limpa' works", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  agg <- prolfqua::AggregateLimpa$new(lfqdata, "protein")
  lfq_agg <- agg$aggregate()

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::build_contrast_analysis(lfq_agg, "~ group_", contrasts, method = "limpa")

  expect_true(inherits(fa, "ContrastsLimpaFacade"))
  res <- fa$get_contrasts()
  expect_equal(nrow(res), dplyr::n_distinct(lfq_agg$data_long()$protein_Id))
  expect_true("diff" %in% colnames(res))
})


test_that("ContrastsLimpaNestedFacade end-to-end (precursor input)", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::ContrastsLimpaNestedFacade$new(lfqdata, "~ group_", contrasts)

  expect_true(inherits(fa, "ContrastsLimpaNestedFacade"))
  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  # one protein-level contrast per protein of the precursor input
  expect_equal(nrow(res), dplyr::n_distinct(lfqdata$data_long()$protein_Id))
  expect_true(all(res$modelName == "limpa_nested"))
  expect_true("diff" %in% colnames(res))
  expect_true("p.value" %in% colnames(res))
  expect_true("FDR" %in% colnames(res))
})


test_that("build_contrast_analysis with method='limpa_nested' works", {
  skip_if_not_installed("limpa")
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 30)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "limpa_nested")

  expect_true(inherits(fa, "ContrastsLimpaNestedFacade"))
  res <- fa$get_contrasts()
  expect_equal(nrow(res), dplyr::n_distinct(lfqdata$data_long()$protein_Id))
})


test_that("ContrastsLimpaFacade accepts aggregated LFQData without opt_se", {
  skip_if_not_installed("limpa")
  # As an "aggregated"-level facade, limpa operates on whatever the normal
  # aggregation pipeline produced. opt_se is optional: when absent (e.g.
  # medpolish-aggregated input), vooma is fit without an SE predictor.
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 20)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq
  agg <- prolfqua::AggregateMedpolish$new(lfqdata, "protein")$aggregate()
  expect_true(length(agg$get_config()$opt_se) == 0 || nchar(agg$get_config()$opt_se) == 0)

  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  fa <- prolfqua::ContrastsLimpaFacade$new(agg, "~ group_", contrasts)
  expect_true(inherits(fa, "ContrastsLimpaFacade"))
  res <- fa$get_contrasts()
  # every medpolish protein is observed in groups A and Ctrl
  expect_equal(nrow(res), dplyr::n_distinct(agg$data_long()$protein_Id))
  expect_true(all(res$modelName == "limpa"))
})

test_that("limpa and limpa_nested equal limpa's own dpcQuant + dpcDE pipeline", {
  skip_if_not_installed("limpa")
  istar <- withr::with_seed(3, sim_lfq_data_peptide_config(n_proteins = 30))
  lfq <- LFQData$new(istar$data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  contrasts <- c(A_vs_Ctrl = "group_A - group_Ctrl")

  # reference: limpa's documented workflow on the same peptide matrix (group A is the baseline level)
  wide <- lfq$data_wide(as.matrix = TRUE)
  protein <- wide$rowdata$protein_Id
  dpc <- limpa::dpc(wide$data)
  protein_elist <- limpa::dpcQuant(wide$data, protein.id = protein, dpc = dpc, dpc.slope = 0.8, verbose = FALSE)
  design <- model.matrix(~group_, data = wide$annotation)
  fit <- limpa::dpcDE(protein_elist, design, plot = FALSE)
  contrast_matrix <- matrix(c(0, 0, -1), ncol = 1, dimnames = list(colnames(design), "A_vs_Ctrl"))
  reference <- limma::topTable(
    limma::eBayes(limma::contrasts.fit(fit, contrast_matrix)),
    coef = 1,
    number = Inf,
    sort.by = "none"
  )

  lfq_agg <- AggregateLimpa$new(lfq, "protein")$aggregate()
  aggregated <- ContrastsLimpaFacade$new(lfq_agg, "~ group_", contrasts)$get_contrasts()
  nested <- ContrastsLimpaNestedFacade$new(lfq, "~ group_", contrasts)$get_contrasts()
  for (ct in list(aggregated, nested)) {
    expect_setequal(ct$protein_Id, rownames(reference))
    ref <- reference[ct$protein_Id, ]
    expect_equal(ct$diff, ref$logFC, tolerance = 1e-8)
    expect_equal(ct$p.value, ref$P.Value, tolerance = 1e-8)
  }
  # a non-default span reaches limpa::voomaLmFitWithImputation()
  fit_span <- limpa::dpcDE(protein_elist, design, plot = FALSE, span = 0.3)
  reference_span <- limma::topTable(
    limma::eBayes(limma::contrasts.fit(fit_span, contrast_matrix)),
    coef = 1,
    number = Inf,
    sort.by = "none"
  )
  with_span <- ContrastsLimpaFacade$new(lfq_agg, "~ group_", contrasts, span = 0.3)$get_contrasts()
  expect_equal(with_span$p.value, reference_span[with_span$protein_Id, "P.Value"], tolerance = 1e-8)
  expect_gt(max(abs(with_span$p.value - reference[with_span$protein_Id, "P.Value"])), 1e-6)
})
