## Tests for ContrastsFacades and build_contrast_analysis

# Shared test data (peptide-level, used for lm / lm_missing / limma / ropeca)
make_peptide_lfqdata <- function(Nprot = 30) {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = Nprot, with_missing = FALSE, seed = 42)

  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq
  keep <- lfqdata$data_long() |>
    dplyr::group_by(protein_Id) |>
    dplyr::summarise(n_peptides = dplyr::n_distinct(peptide_Id), .groups = "drop") |>
    dplyr::filter(n_peptides > 1)
  lfqdata$set_data(dplyr::semi_join(lfqdata$data_long(), keep, by = "protein_Id"))
  lfqdata
}

# Protein-level data with nr_peptides (used for deqms)
make_protein_lfqdata <- function(Nprot = 30) {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot, seed = 42)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  list(lfqdata = lfqdata, raw = istar)
}

CONTRASTS <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
MODELSTR <- "~ group_"
MODELSTR_LMER <- "~ group_ + (1 | peptide_Id) + (1 | sampleName)"

# Helper: check that a facade has the required interface
check_facade_interface <- function(fa) {
  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  expected_cols <- c("modelName", "estimate_type", "contrast", "diff", "FDR", "p.value")
  for (col in expected_cols) {
    expect_true(col %in% colnames(res), info = paste("missing column:", col))
  }
  # clean break: the legacy `facade` column is gone; modelName carries identity
  expect_false("facade" %in% colnames(res))
  # get_contrasts must not contain NA diff
  expect_true(!any(is.na(res$diff)), info = "get_contrasts() must not return NA diff rows")

  missing <- fa$get_missing()
  expect_true(is.data.frame(missing))
  expect_true("contrast" %in% colnames(missing))

  wide <- fa$to_wide()
  expect_true(is.data.frame(wide))
  expect_true(nrow(wide) > 0)

  plotter <- fa$get_Plotter()
  expect_true(!is.null(plotter))
}

# Regression: get_contrasts(), to_wide(), and get_Plotter() must all report the
# selected facade key, not the inner model name (e.g. "WaldTest"). lm wraps
# Contrasts -> ContrastsModerated, the same path that previously leaked.
test_that("facade modelName is consistent across get_contrasts/to_wide/get_Plotter", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "lm")

  expect_setequal(unique(fa$get_contrasts()$modelName), "lm")

  wide <- fa$to_wide(columns = "modelName")
  wide_modelname <- unique(unlist(wide[grep("^modelName", names(wide))]))
  expect_setequal(wide_modelname, "lm")

  pl <- fa$get_Plotter()
  pl_df <- pl$.__enclos_env__$private$.contrast_df
  expect_setequal(unique(pl_df$modelName), "lm")
})

# ---- ContrastsLimmaFacade ----

test_that("ContrastsLimmaFacade initialises and returns correct structure", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::ContrastsLimmaFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsLimmaFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  check_facade_interface(fa)
  expect_true(all(fa$get_contrasts()$modelName == "limma"))
})

# ---- ContrastsLMFacade ----

test_that("ContrastsLMFacade initialises and returns correct structure", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::ContrastsLMFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsLMFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  check_facade_interface(fa)
  expect_true(all(fa$get_contrasts()$modelName == "lm"))
})

# ---- ContrastsRLMFacade ----

test_that("ContrastsRLMFacade initialises and returns correct structure", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::ContrastsRLMFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsRLMFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  check_facade_interface(fa)
  expect_true(all(fa$get_contrasts()$modelName == "rlm"))
})

# ---- ContrastsLmerNestedFacade ----

test_that("ContrastsLmerNestedFacade initialises and returns correct structure", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::ContrastsLmerNestedFacade$new(lfqdata, MODELSTR_LMER, CONTRASTS)

  expect_true(inherits(fa, "ContrastsLmerNestedFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  check_facade_interface(fa)
  expect_true(all(fa$get_contrasts()$modelName == "lmer_nested"))
})

test_that("ContrastsLmerNestedFacade augments fixed-effects-only modelstr with random effects", {
  # Caller passes only the fixed effects "~ group_"; facade derives
  # (1 | <deepest child key>) + (1 | <sample name>) from the LFQData.
  lfqdata <- make_peptide_lfqdata()
  fa_auto <- prolfqua::ContrastsLmerNestedFacade$new(lfqdata, "~ group_", CONTRASTS)
  fa_full <- prolfqua::ContrastsLmerNestedFacade$new(lfqdata, MODELSTR_LMER, CONTRASTS)

  expect_true(inherits(fa_auto$model$model_strategy, "StrategyLmer"))
  # Both formulas should be identical after augmentation.
  expect_equal(
    deparse1(fa_auto$model$model_strategy$formula),
    deparse1(fa_full$model$model_strategy$formula)
  )
  # And the contrast results should match.
  expect_equal(
    fa_auto$get_contrasts()$diff,
    fa_full$get_contrasts()$diff
  )
})

# ---- ContrastsLMMissingFacade ----

test_that("ContrastsLMMissingFacade initialises and returns correct structure", {
  # ContrastsMissing requires protein-level data (hierarchyDepth == len(hierarchy_keys()))
  obj <- make_protein_lfqdata()
  lfqdata <- obj$lfqdata
  fa <- prolfqua::ContrastsLMMissingFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsLMMissingFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  expect_true(!is.null(fa$missing_contrast))
  expect_true(!is.null(fa$merged))
  check_facade_interface(fa)
  expect_true(all(fa$get_contrasts()$modelName == "lm_missing"))
})

# ---- ContrastsDEqMSFacade ----

test_that("ContrastsDEqMSFacade initialises and returns correct structure", {
  obj <- make_protein_lfqdata()
  lfqdata <- obj$lfqdata
  fa <- prolfqua::ContrastsDEqMSFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsDEqMSFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  expect_equal(fa$contrast$count_column, "nr_peptides")
  check_facade_interface(fa)
  expect_true(all(fa$get_contrasts()$modelName == "deqms"))
})

# ---- ContrastsROPECANestedFacade ----

test_that("ContrastsROPECANestedFacade initialises and returns correct structure", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::ContrastsROPECANestedFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsROPECANestedFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))

  # Facade normalises ROPECA output to standard column names
  check_facade_interface(fa)
  expect_true(all(fa$get_contrasts()$modelName == "ropeca_nested"))

  res <- fa$get_contrasts()
  # Heuristically derived columns should have real values (not all NA)
  expect_true(!all(is.na(res$std.error)))
  expect_true(!all(is.na(res$df)))
  expect_true(!all(is.na(res$sigma)))
  expect_true(!all(is.na(res$conf.low)))
  expect_true(!all(is.na(res$conf.high)))
  # std.error = diff / statistic
  nonzero <- res$statistic != 0 & !is.na(res$statistic)
  expect_equal(res$std.error[nonzero], res$diff[nonzero] / res$statistic[nonzero])
  # df = number of contributing peptides (positive integers)
  expect_true(all(res$df > 0 & res$df == floor(res$df)))
  # ROPECA raw columns should NOT appear in facade output
  expect_false("beta.based.significance" %in% colnames(res))
  expect_false("FDR.beta.based.significance" %in% colnames(res))
  expect_false("mad.estimate" %in% colnames(res))
  expect_false("n_not_na" %in% colnames(res))
})

# ---- build_contrast_analysis ----

test_that("build_contrast_analysis dispatches to ContrastsLMFacade for method='lm'", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "lm")
  expect_true(inherits(fa, "ContrastsLMFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsRLMFacade for method='rlm'", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "rlm")
  expect_true(inherits(fa, "ContrastsRLMFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsLmerNestedFacade for method='lmer_nested'", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR_LMER, CONTRASTS, method = "lmer_nested")
  expect_true(inherits(fa, "ContrastsLmerNestedFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsLMMissingFacade for method='lm_missing'", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "lm_missing")
  expect_true(inherits(fa, "ContrastsLMMissingFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsLimmaFacade for method='limma'", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "limma")
  expect_true(inherits(fa, "ContrastsLimmaFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsDEqMSFacade for method='deqms'", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "deqms")
  expect_true(inherits(fa, "ContrastsDEqMSFacade"))
  expect_equal(fa$contrast$count_column, "nr_peptides")
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsROPECANestedFacade for method='ropeca_nested'", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "ropeca_nested")
  expect_true(inherits(fa, "ContrastsROPECANestedFacade"))
  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  # Standard column names, not ROPECA-specific ones
  expect_true("p.value" %in% colnames(res))
  expect_true("FDR" %in% colnames(res))
  expect_false("beta.based.significance" %in% colnames(res))
})

test_that("build_contrast_analysis dispatches to ContrastsFirthNestedFacade for method='firth_nested'", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "firth_nested")
  expect_true(inherits(fa, "ContrastsFirthNestedFacade"))
  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(all(res$modelName == "firth_nested"))
})

test_that("Firth nested preparation completes sparse LFQData before encoding missingness", {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 20, weight_missing = 1, seed = 2)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  response <- lfqdata$response()

  lfq_sparse <- lfqdata$get_copy()
  lfq_sparse$set_data(lfq_sparse$data_long()[!is.na(lfq_sparse$data_long()[[response]]), ])
  expect_setequal(unique(prolfqua::encode_bin_resp(lfq_sparse)$bin_resp), 1L)

  lfq_missing <- prolfqua:::.prepare_logistf_lfqdata(lfq_sparse)
  expect_setequal(unique(lfq_missing$data_long()$bin_resp), c(0L, 1L))
})

test_that("build_contrast_analysis defaults to lm when no method specified", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS)
  expect_true(inherits(fa, "ContrastsLMFacade"))
})

test_that("aggregated facades error on peptide-level LFQData", {
  lfqdata <- make_peptide_lfqdata()

  expect_error(
    prolfqua::ContrastsLMFacade$new(lfqdata, MODELSTR, CONTRASTS),
    "requires aggregated LFQData"
  )
  expect_error(
    prolfqua::ContrastsRLMFacade$new(lfqdata, MODELSTR, CONTRASTS),
    "requires aggregated LFQData"
  )
  expect_error(
    prolfqua::ContrastsLimmaFacade$new(lfqdata, MODELSTR, CONTRASTS),
    "requires aggregated LFQData"
  )
  expect_error(
    prolfqua::ContrastsLMMissingFacade$new(lfqdata, MODELSTR, CONTRASTS),
    "requires aggregated LFQData"
  )
  expect_error(
    prolfqua::ContrastsDEqMSFacade$new(lfqdata, MODELSTR, CONTRASTS),
    "requires aggregated LFQData"
  )
})

# ---- get_missing with data that has missingness ----

test_that("get_missing returns non-empty for data with missing groups", {
  # Simulate with missingness so some proteins can't be estimated
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 80, seed = 42)

  lfq_peptide <- prolfqua::LFQData$new(istar$data, istar$config)
  lfq_peptide <- lfq_peptide$get_Transformer()$log2()$lfq
  aggregator <- prolfqua::AggregateMedpolish$new(lfq_peptide, "protein")
  lfq_protein <- aggregator$aggregate()

  contrasts2 <- c("A_vs_Ctrl" = "group_A - group_Ctrl", "B_vs_Ctrl" = "group_B - group_Ctrl")

  fa_lm <- prolfqua::ContrastsLMFacade$new(lfq_protein, MODELSTR, contrasts2)
  fa_limma <- prolfqua::ContrastsLimmaFacade$new(lfq_protein, MODELSTR, contrasts2)
  fa_lm_missing <- prolfqua::ContrastsLMMissingFacade$new(lfq_protein, MODELSTR, contrasts2)

  # lm should have missing proteins (it drops unfittable ones)
  lm_missing <- fa_lm$get_missing()
  expect_true(is.data.frame(lm_missing))
  expect_true("protein_Id" %in% colnames(lm_missing))
  expect_true("contrast" %in% colnames(lm_missing))

  # limma should also have missing proteins (NA diff rows are filtered)
  limma_missing <- fa_limma$get_missing()
  expect_true(is.data.frame(limma_missing))
  expect_true("protein_Id" %in% colnames(limma_missing))

  # limma get_contrasts must have no NA diff
  limma_res <- fa_limma$get_contrasts()
  expect_true(!any(is.na(limma_res$diff)))

  # lm_missing should have fewer (or no) missing proteins thanks to imputation
  lm_missing_missing <- fa_lm_missing$get_missing()
  expect_true(nrow(lm_missing_missing) <= nrow(lm_missing))

  # get_missing + get_contrasts should cover all input proteins × contrasts
  n_input_proteins <- dplyr::n_distinct(lfq_protein$data_long()$protein_Id)
  n_contrasts <- length(contrasts2)
  for (fa in list(fa_lm, fa_limma, fa_lm_missing)) {
    n_estimated <- nrow(fa$get_contrasts())
    n_missing <- nrow(fa$get_missing())
    expect_equal(n_estimated + n_missing, n_input_proteins * n_contrasts)
  }
})

test_that("get_missing returns empty data.frame when all proteins are estimable", {
  # No missingness — all proteins should be estimable
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::ContrastsLMFacade$new(lfqdata, MODELSTR, CONTRASTS)
  missing <- fa$get_missing()
  expect_true(is.data.frame(missing))
  expect_equal(nrow(missing), 0)
})

test_that("nested facades error on aggregated LFQData", {
  lfqdata <- make_protein_lfqdata()$lfqdata

  expect_error(
    prolfqua::ContrastsROPECANestedFacade$new(lfqdata, MODELSTR, CONTRASTS),
    "requires LFQData with additional hierarchy below `subject_id\\(\\)`"
  )
  expect_error(
    prolfqua::ContrastsLmerNestedFacade$new(lfqdata, MODELSTR_LMER, CONTRASTS),
    "requires LFQData with additional hierarchy below `subject_id\\(\\)`"
  )
})
