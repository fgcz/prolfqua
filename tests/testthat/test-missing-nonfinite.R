# Missing and non-finite values: NA intensities are missing observations; Inf, -Inf, NaN responses, NA keys and
# empty inputs raise classed errors before fitting or aggregation; log transforms turn 0 (not detected) into NA.

md_peptide <- function(with_missing = TRUE) {
  sim <- sim_lfq_data_peptide_config(Nprot = 10, with_missing = with_missing, seed = 42)
  LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
}

md_protein <- function(with_missing = TRUE) {
  md_peptide(with_missing)$get_Aggregator()$aggregate()
}

md_replace <- function(lfq, column, value, rows = 1L) {
  data <- lfq$data_long()
  data[[column]][rows] <- value
  LFQData$new(data, lfq$get_config())
}

MD_CONTRASTS <- c(A_vs_Ctrl = "group_A - group_Ctrl")

test_that("Inf, -Inf and NaN responses raise before fitting and aggregation", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- md_protein()
  peptide <- md_peptide()
  for (value in c(Inf, -Inf, NaN)) {
    bad_protein <- md_replace(protein, protein$response(), value)
    expect_error(
      build_contrast_analysis(bad_protein, "~ group_", MD_CONTRASTS, "lm"),
      class = "prolfqua_error_nonfinite_response"
    )
    expect_error(
      build_model(bad_protein, strategy_lm(paste(bad_protein$response(), "~ group_"))),
      class = "prolfqua_error_nonfinite_response"
    )
    bad_peptide <- md_replace(peptide, peptide$response(), value)
    for (method in c("medpolish", "rlm", "topN")) {
      expect_error(bad_peptide$get_Aggregator(method), class = "prolfqua_error_nonfinite_response")
    }
  }
  expect_error(
    build_contrast_analysis(md_replace(protein, protein$response(), Inf), "~ group_", MD_CONTRASTS, "lm"),
    regexp = protein$data_long()$protein_Id[1],
    fixed = TRUE
  )
})

test_that("NA in hierarchy, factor or sample keys raises", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- md_protein()
  for (key in c("protein_Id", "group_", protein$sample_name())) {
    expect_error(
      build_contrast_analysis(md_replace(protein, key, NA), "~ group_", MD_CONTRASTS, "lm"),
      class = "prolfqua_error_missing_key",
      regexp = key
    )
  }
})

test_that("an input without rows raises", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- md_protein()
  empty <- LFQData$new(protein$data_long()[0, ], protein$get_config())
  expect_error(
    build_contrast_analysis(empty, "~ group_", MD_CONTRASTS, "lm"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("a zero intensity becomes NA after the log transform and keeps its protein", {
  sim <- sim_lfq_data_peptide_config(Nprot = 10, with_missing = FALSE, seed = 42)
  sim$data[[sim$config$get_response()]][1] <- 0
  expect_message(
    peptide <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq,
    "1 zeros.*set to NA"
  )
  expect_true(is.na(peptide$data_long()[[peptide$response()]][1]))
  withr::local_options(prolfqua.verbose = FALSE)
  contrasts <- build_contrast_analysis(
    peptide$get_Aggregator()$aggregate(),
    "~ group_",
    MD_CONTRASTS,
    "lm"
  )$get_contrasts()
  expect_true(sim$data$protein_Id[1] %in% contrasts$protein_Id)
})

test_that("a protein with an all-NA response is absent from the contrasts and listed as missing", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- md_protein(with_missing = FALSE)
  target <- unique(protein$data_long()$protein_Id)[1]
  all_na <- md_replace(protein, protein$response(), NA_real_, protein$data_long()$protein_Id == target)
  for (method in c("lm", "limma")) {
    facade <- muffle_warnings_matching(
      build_contrast_analysis(all_na, "~ group_", MD_CONTRASTS, method),
      paste(one_group_warning, partial_na_warning, sep = "|")
    )
    expect_false(target %in% facade$get_contrasts()$protein_Id, label = method)
    expect_true(target %in% facade$get_missing()$protein_Id, label = method)
  }
})

test_that("a protein with a constant response gives finite or NA results, no error", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- md_protein(with_missing = FALSE)
  target <- unique(protein$data_long()$protein_Id)[1]
  constant <- md_replace(protein, protein$response(), 20, protein$data_long()$protein_Id == target)
  for (method in c("lm", "limma")) {
    contrasts <- muffle_warnings_matching(
      build_contrast_analysis(constant, "~ group_", MD_CONTRASTS, method)$get_contrasts(),
      paste("essentially perfect fit", partial_na_warning, sep = "|")
    )
    row <- contrasts[contrasts$protein_Id == target, c("diff", "std.error", "p.value")]
    expect_equal(nrow(row), 1L, label = method)
    expect_false(any(is.infinite(unlist(row)) | is.nan(unlist(row))), label = method)
  }
})

#' @srrstats {G5.3} Contrast tables of the intensity methods lm, rlm, limma, lm_impute, limma_impute and deqms on data
#'   with missing values contain no Inf or NaN; NA marks a contrast that cannot be estimated. Detection models may
#'   return NaN p-values for proteins detected in every sample (Modelling methods vignette).
#' @noRd
NULL

test_that("contrast tables contain no Inf or NaN", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- md_protein(with_missing = TRUE)
  for (method in c("lm", "rlm", "limma", "lm_impute", "limma_impute", "deqms")) {
    contrasts <- muffle_warnings_matching(
      build_contrast_analysis(protein, "~ group_", MD_CONTRASTS, method)$get_contrasts(),
      paste(
        one_group_warning,
        partial_na_warning,
        loess_ties_warning,
        rlm_convergence_warning,
        "essentially perfect fit",
        sep = "|"
      )
    )
    numeric <- Filter(is.numeric, contrasts)
    expect_gt(nrow(contrasts), 0L)
    expect_false(any(vapply(numeric, function(x) any(is.infinite(x) | is.nan(x)), logical(1))), label = method)
  }
})

test_that("moderation with equal variances ignores an NA residual df for the prior df", {
  contrast_df <- data.frame(
    sigma = rep(0.3, 5),
    df = c(6, 6, NA, 6, 6),
    statistic = c(2, -1, 0.5, 3, -2),
    std.error = rep(0.3, 5),
    diff = c(0.6, -0.3, 0.15, 0.9, -0.6)
  )
  res <- moderated_p_limma(contrast_df)
  expect_true(all(is.finite(res$moderated.df.prior)))
  expect_true(all(is.finite(res$moderated.p.value[!is.na(contrast_df$df)])))
})

test_that("which_missing rejects NA input", {
  expect_error(which_missing(c(1, NA, 3)), class = "prolfqua_error_bad_argument")
})

test_that("simulators with with_missing = FALSE give complete data with positive peptide counts", {
  sims <- list(
    peptide = function(m) sim_lfq_data_peptide_config(Nprot = 10, with_missing = m, seed = 42),
    protein = function(m) sim_lfq_data_protein_config(Nprot = 10, with_missing = m, seed = 42),
    two_factor = function(m) sim_lfq_data_2factor_config(Nprot = 10, with_missing = m, seed = 42)
  )
  for (name in names(sims)) {
    complete <- sims[[name]](FALSE)
    with_na <- sims[[name]](TRUE)
    counts <- complete$data[[complete$config$nr_children]]
    expect_true(all(counts > 0), label = name)
    expect_false(anyNA(complete$data$abundance), label = name)
    observed <- with_na$data[!is.na(with_na$data$abundance), ]
    expect_lt(nrow(observed), nrow(complete$data))
    key <- c(complete$config$sample_name, complete$config$hierarchy_keys())
    joined <- dplyr::inner_join(observed, complete$data, by = key, suffix = c(".na", ".complete"))
    expect_equal(nrow(joined), nrow(observed), label = name)
    expect_equal(joined$abundance.na, joined$abundance.complete, label = name)
  }
})
