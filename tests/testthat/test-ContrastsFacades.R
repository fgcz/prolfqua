## Tests for ContrastsFacades and build_contrast_analysis

# Shared test data (peptide-level, used for lm / lm_missing / limma / ropeca)
make_peptide_lfqdata <- function(Nprot = 30) {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = Nprot, seed = 42)
  istar$config <- prolfqua::old2new(istar$config)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$get_Transformer()$log2()$lfq
}

# Protein-level data with nr_peptides (used for deqms)
make_protein_lfqdata <- function(Nprot = 30) {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot, seed = 42)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  list(lfqdata = lfqdata, raw = istar)
}

CONTRASTS <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
MODELSTR  <- "~ group_"

# Helper: check that a facade has the required interface
check_facade_interface <- function(fa) {
  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  expected_cols <- c("modelName", "contrast", "diff", "FDR", "p.value")
  for (col in expected_cols) {
    expect_true(col %in% colnames(res), info = paste("missing column:", col))
  }

  wide <- fa$to_wide()
  expect_true(is.data.frame(wide))
  expect_true(nrow(wide) > 0)

  plotter <- fa$get_Plotter()
  expect_true(!is.null(plotter))
}

# ---- ContrastsLimmaFacade ----

test_that("ContrastsLimmaFacade initialises and returns correct structure", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::ContrastsLimmaFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsLimmaFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  check_facade_interface(fa)
})

# ---- ContrastsLMFacade ----

test_that("ContrastsLMFacade initialises and returns correct structure", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::ContrastsLMFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsLMFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  check_facade_interface(fa)
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
})

# ---- ContrastsDEqMSFacade ----

test_that("ContrastsDEqMSFacade initialises and returns correct structure", {
  obj <- make_protein_lfqdata()
  lfqdata <- obj$lfqdata
  istar   <- obj$raw
  count_df <- dplyr::select(istar$data,
    dplyr::all_of(c(istar$config$hierarchy_keys_depth(), "nr_peptides"))) |>
    dplyr::distinct()

  fa <- prolfqua::ContrastsDEqMSFacade$new(lfqdata, MODELSTR, CONTRASTS,
                                            count_df = count_df,
                                            count_column = "nr_peptides")

  expect_true(inherits(fa, "ContrastsDEqMSFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))
  check_facade_interface(fa)
})

# ---- ContrastsROPECAFacade ----

test_that("ContrastsROPECAFacade initialises and returns correct structure", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::ContrastsROPECAFacade$new(lfqdata, MODELSTR, CONTRASTS)

  expect_true(inherits(fa, "ContrastsROPECAFacade"))
  expect_true(!is.null(fa$model))
  expect_true(!is.null(fa$contrast))

  # Facade normalises ROPECA output to standard column names
  check_facade_interface(fa)

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
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "lm")
  expect_true(inherits(fa, "ContrastsLMFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsLMMissingFacade for method='lm_missing'", {
  lfqdata <- make_protein_lfqdata()$lfqdata
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "lm_missing")
  expect_true(inherits(fa, "ContrastsLMMissingFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsLimmaFacade for method='limma'", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "limma")
  expect_true(inherits(fa, "ContrastsLimmaFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsDEqMSFacade for method='deqms'", {
  obj <- make_protein_lfqdata()
  lfqdata <- obj$lfqdata
  istar   <- obj$raw
  count_df <- dplyr::select(istar$data,
    dplyr::all_of(c(istar$config$hierarchy_keys_depth(), "nr_peptides"))) |>
    dplyr::distinct()

  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS,
                                           method = "deqms",
                                           count_df = count_df,
                                           count_column = "nr_peptides")
  expect_true(inherits(fa, "ContrastsDEqMSFacade"))
  check_facade_interface(fa)
})

test_that("build_contrast_analysis dispatches to ContrastsROPECAFacade for method='ropeca'", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "ropeca")
  expect_true(inherits(fa, "ContrastsROPECAFacade"))
  res <- fa$get_contrasts()
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  # Standard column names, not ROPECA-specific ones
  expect_true("p.value" %in% colnames(res))
  expect_true("FDR" %in% colnames(res))
  expect_false("beta.based.significance" %in% colnames(res))
})

test_that("build_contrast_analysis errors when deqms missing count_df", {
  lfqdata <- make_peptide_lfqdata()
  expect_error(
    prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS, method = "deqms"),
    "count_df and count_column are required"
  )
})

test_that("build_contrast_analysis defaults to lm when no method specified", {
  lfqdata <- make_peptide_lfqdata()
  fa <- prolfqua::build_contrast_analysis(lfqdata, MODELSTR, CONTRASTS)
  expect_true(inherits(fa, "ContrastsLMFacade"))
})
