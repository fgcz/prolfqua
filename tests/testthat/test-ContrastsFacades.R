## Tests for ContrastsFacades and build_contrast_analysis

# Shared test data (peptide-level, used for lm / limma / ropeca)
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

# Protein-level data with nrPeptides (used for deqms)
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

# ---- ContrastsDEqMSFacade ----

test_that("ContrastsDEqMSFacade initialises and returns correct structure", {
  obj <- make_protein_lfqdata()
  lfqdata <- obj$lfqdata
  fa <- prolfqua::ContrastsDEqMSFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsDEqMSFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  expect_equal(fa$contrast$count_column, "nrPeptides")
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
  expect_equal(fa$contrast$count_column, "nrPeptides")
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

  lfq_missing <- prolfqua:::.prepare_detection_lfqdata(lfq_sparse)
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

  # lm should have missing proteins (it drops unfittable ones)
  missing_lm <- fa_lm$get_missing()
  expect_true(is.data.frame(missing_lm))
  expect_true("protein_Id" %in% colnames(missing_lm))
  expect_true("contrast" %in% colnames(missing_lm))

  # limma should also have missing proteins (NA diff rows are filtered)
  limma_missing <- fa_limma$get_missing()
  expect_true(is.data.frame(limma_missing))
  expect_true("protein_Id" %in% colnames(limma_missing))

  # limma get_contrasts must have no NA diff
  limma_res <- fa_limma$get_contrasts()
  expect_true(!any(is.na(limma_res$diff)))

  # get_missing + get_contrasts should cover all input proteins × contrasts
  n_input_proteins <- dplyr::n_distinct(lfq_protein$data_long()$protein_Id)
  n_contrasts <- length(contrasts2)
  for (fa in list(fa_lm, fa_limma)) {
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

test_that("ContrastsROPECANestedFacade get_contrasts(all = TRUE) appends the ROPECA columns", {
  fa <- prolfqua::ContrastsROPECANestedFacade$new(make_peptide_lfqdata(), MODELSTR, CONTRASTS)
  standard <- colnames(fa$get_contrasts())
  everything <- colnames(fa$get_contrasts(all = TRUE))
  expect_equal(everything[seq_along(standard)], standard)
  expect_gt(length(everything), length(standard))
})

test_that("ContrastsTable$new stores a passed ContrastConfiguration", {
  fa <- prolfqua::ContrastsLMFacade$new(make_protein_lfqdata()$lfqdata, MODELSTR, CONTRASTS)
  cfg <- prolfqua::ContrastConfiguration$new(subject_id = "protein_Id", fdr_col = "FDR")
  tab <- prolfqua::ContrastsTable$new(fa$get_contrasts(), subject_id = "protein_Id", config = cfg)
  expect_identical(tab$get_config(), cfg)
  expect_error(
    prolfqua::ContrastsTable$new(fa$get_contrasts(), config = list()),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("a contrast naming a group without samples is a typed error", {
  lfq <- make_protein_lfqdata()$lfqdata
  lfq$set_data(dplyr::filter(lfq$data_long(), .data$group_ != "A"))
  expect_error(
    prolfqua::ContrastsLMFacade$new(lfq, MODELSTR, c(AvsCtrl = "group_A - group_Ctrl")),
    regexp = "group_A",
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    prolfqua::build_contrast_analysis(lfq, MODELSTR, c(AvsCtrl = "group_A - group_Ctrl"), method = "limma"),
    class = "prolfqua_error_bad_argument"
  )
  # contrasts between groups with data still work
  fa <- prolfqua::ContrastsLMFacade$new(lfq, MODELSTR, c(BvsCtrl = "group_B - group_Ctrl"))
  expect_gt(nrow(fa$get_contrasts()), 0)
})

test_that("the firth facade keeps a protein key not named protein_Id", {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 20, weight_missing = 0.5, seed = 42)
  data <- dplyr::rename(istar$data, prot = "protein_Id")
  config <- istar$config$clone(deep = TRUE)
  names(config$hierarchy) <- "prot"
  lfq <- prolfqua::LFQData$new(data, config)
  fa <- prolfqua::ContrastsFirthFacade$new(lfq, MODELSTR, CONTRASTS)
  expect_equal(fa$contrast$subject_id, "prot")
  res <- fa$get_contrasts()
  expect_true("prot" %in% colnames(res))
  expect_false("protein_Id" %in% colnames(res))
})

test_that(".compute_missing keys missing pairs by the configured contrast column", {
  lfq <- make_protein_lfqdata()$lfqdata
  proteins <- unique(lfq$data_long()$protein_Id)
  result <- data.frame(protein_Id = proteins[-1], Bait = "AvsCtrl")
  missing <- prolfqua:::.compute_missing(lfq, "AvsCtrl", result, contrast_col = "Bait")
  expect_equal(colnames(missing), c("protein_Id", "Bait"))
  expect_equal(missing$protein_Id, proteins[1])
})
