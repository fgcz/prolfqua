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

# Proteins with zero residual variance: one constant (all 20), one noiseless with a real difference (A 22, rest 20)
md_degenerate <- function() {
  protein <- md_protein(with_missing = FALSE)
  data <- protein$data_long()
  response <- protein$response()
  ids <- unique(data$protein_Id)
  data[[response]][data$protein_Id == ids[1]] <- 20
  noiseless <- data$protein_Id == ids[2]
  data[[response]][noiseless] <- ifelse(data$group_[noiseless] == "A", 22, 20)
  list(lfq = LFQData$new(data, protein$get_config()), constant = ids[1], noiseless = ids[2], ids = ids)
}

md_contrasts_counting_zero_variance <- function(lfq, method) {
  n_zero <- 0L
  contrasts <- withCallingHandlers(
    muffle_warnings_matching(
      build_contrast_analysis(lfq, "~ group_", MD_CONTRASTS, method)$get_contrasts(),
      paste("essentially perfect fit", loess_ties_warning, rlm_convergence_warning, sep = "|")
    ),
    prolfqua_zero_variance = function(w) {
      n_zero <<- n_zero + 1L
      invokeRestart("muffleWarning")
    }
  )
  list(contrasts = contrasts, n_zero = n_zero)
}

test_that("zero residual variance keeps the difference and sets the test to NA, with one warning", {
  withr::local_options(prolfqua.verbose = FALSE)
  degenerate <- md_degenerate()
  for (method in c("lm", "rlm", "limma", "deqms", "lm_impute", "limma_impute")) {
    res <- md_contrasts_counting_zero_variance(degenerate$lfq, method)
    expect_gte(res$n_zero, 1L)
    rows <- res$contrasts[res$contrasts$protein_Id %in% c(degenerate$constant, degenerate$noiseless), ]
    expect_true(all(is.na(rows$p.value)), label = method)
    expect_true(all(is.na(rows$std.error)), label = method)
    constant <- rows[rows$protein_Id == degenerate$constant, ]
    expect_equal(constant$diff, 0, tolerance = 1e-8, label = method)
    if (method != "rlm") {
      noiseless <- rows[rows$protein_Id == degenerate$noiseless, ]
      expect_equal(noiseless$diff, 2, tolerance = 1e-8, label = method)
    }
  }
})

test_that("proteins with zero residual variance do not change the variance prior of the others", {
  withr::local_options(prolfqua.verbose = FALSE)
  degenerate <- md_degenerate()
  others <- setdiff(degenerate$ids, c(degenerate$constant, degenerate$noiseless))
  without <- LFQData$new(
    degenerate$lfq$data_long()[degenerate$lfq$data_long()$protein_Id %in% others, ],
    degenerate$lfq$get_config()
  )
  for (method in c("lm", "limma")) {
    with_degenerate <- md_contrasts_counting_zero_variance(degenerate$lfq, method)$contrasts
    reference <- md_contrasts_counting_zero_variance(without, method)$contrasts
    with_degenerate <- with_degenerate[match(reference$protein_Id, with_degenerate$protein_Id), ]
    expect_equal(with_degenerate$p.value, reference$p.value, tolerance = 1e-10, label = method)
  }
})

test_that("every subject with zero residual variance gives NA tests, no error", {
  withr::local_options(prolfqua.verbose = FALSE)
  protein <- md_protein(with_missing = FALSE)
  data <- protein$data_long()
  data[[protein$response()]] <- ifelse(data$group_ == "A", 22, 20)
  noiseless <- LFQData$new(data, protein$get_config())
  for (method in c("lm", "limma", "deqms")) {
    res <- md_contrasts_counting_zero_variance(noiseless, method)
    expect_equal(res$n_zero, 1L, label = method)
    expect_true(all(is.na(res$contrasts$p.value)), label = method)
    expect_equal(unname(res$contrasts$diff), rep(2, nrow(res$contrasts)), tolerance = 1e-8, label = method)
  }
})

test_that("moderated_p_limma leaves NA variances out of the prior", {
  contrast_df <- data.frame(
    sigma = c(0.25, 0.32, 0.28, 0.40, 0.35, NA),
    df = rep(6, 6),
    statistic = c(2.1, -1.8, 0.5, 3.0, -2.2, NA),
    std.error = c(0.38, 0.33, 0.4, 0.4, 0.41, NA),
    diff = c(0.8, -0.6, 0.2, 1.2, -0.9, 2)
  )
  with_na <- moderated_p_limma(contrast_df)
  reference <- moderated_p_limma(contrast_df[1:5, ])
  expect_equal(with_na$moderated.p.value[1:5], reference$moderated.p.value)
  expect_true(is.na(with_na$moderated.p.value[6]))
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
