# Test Contrasts R6 classes

test_that("Contrasts (Wald test)", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  muffle_warnings_matching(single_peptide_warning, expr = {
    istar <- sim_lfq_data_peptide_config(n_proteins = 20)
    config <- istar$config

    modelFunction <-
      strategy_lmer("abundance ~ group_ + (1 | peptide_Id) + (1 | sample)")

    mod <- build_model(
      istar$data,
      modelFunction,
      subject_id = config$hierarchy_keys_depth()
    )

    Contr <- c(
      "groupA_vs_Ctrl" = "group_A - group_Ctrl",
      "groupB_vs_Ctrl" = "group_B - group_Ctrl"
    )

    contrastX <- Contrasts$new(mod, Contr)

    # get_linfct
    lf <- contrastX$get_linfct(avg = FALSE)
    expect_true(is.matrix(lf) || is.data.frame(lf))

    # get_contrasts returns data.frame with expected columns
    x <- contrastX$get_contrasts()
    expect_s3_class(x, "data.frame")
    expect_true(all(c("contrast", "diff", "statistic", "p.value", "FDR") %in% colnames(x)))
    expect_equal(x$std.error.unmoderated, x$std.error)
    expect_equal(x$df.unmoderated, x$df)
    expect_true(all(x$p.value > 0 & x$p.value < 1))
    expect_true(all(x$FDR > 0 & x$FDR <= 1))

    # get_contrast_sides
    cs <- contrastX$get_contrast_sides()
    expect_s3_class(cs, "data.frame")
    expect_true(all(c("group_1", "group_2") %in% colnames(cs)))

    # column_description
    cd <- contrastX$column_description()
    expect_s3_class(cd, "data.frame")

    # to_wide
    tw <- contrastX$to_wide()
    expect_s3_class(tw, "data.frame")

    # get_Plotter
    pl <- contrastX$get_Plotter()
    expect_s3_class(pl, "ContrastsPlotter")
  })
})

test_that("ContrastsModerated", {
  istar <- sim_lfq_data_protein_config(n_proteins = 20)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  modelFunction <- strategy_lm("transformedIntensity ~ group_")
  mod <- build_model(lProt, modelFunction)

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  contrast <- Contrasts$new(mod, Contr)
  moderated <- ContrastsModerated$new(contrast)

  raw <- contrast$get_contrasts()
  x <- moderated$get_contrasts()
  expect_s3_class(x, "data.frame")
  expect_true(all(c("diff", "p.value", "FDR") %in% colnames(x)))
  expect_equal(x$std.error.unmoderated, raw$std.error)
  expect_equal(x$df.unmoderated, raw$df)
  expect_equal(x$statistic, x$diff / x$std.error, tolerance = 1e-12)
  expect_equal(
    x$p.value,
    2 * stats::pt(abs(x$statistic), df = x$df, lower.tail = FALSE),
    tolerance = 1e-12
  )
  expect_equal(
    x$conf.high - x$diff,
    -stats::qt(0.025, df = x$df) * x$std.error,
    tolerance = 1e-12
  )

  tw <- moderated$to_wide()
  expect_s3_class(tw, "data.frame")

  pl <- moderated$get_Plotter()
  expect_s3_class(pl, "ContrastsPlotter")
})

test_that("ContrastsROPECA", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 20)
  config <- istar$config$clone(deep = TRUE)
  config$hierarchy_depth <- 2

  modelFunction <- strategy_lm("abundance ~ group_")
  mod <- build_model(
    istar$data,
    modelFunction,
    subject_id = config$hierarchy_keys_depth()
  )

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  contr <- Contrasts$new(mod, Contr)
  ropeca <- ContrastsROPECA$new(contr)

  x <- ropeca$get_contrasts()
  expect_s3_class(x, "data.frame")

  tw <- ropeca$to_wide()
  expect_s3_class(tw, "data.frame")

  pl <- ropeca$get_Plotter()
  expect_s3_class(pl, "ContrastsPlotter")

  # Delegates
  expect_s3_class(ropeca$get_contrast_sides(), "data.frame")
})

test_that("summary_ropeca_median_p_scaled counts only peptides with a p-value", {
  peptides <- data.frame(
    contrast = "c1",
    protein_Id = "P1",
    diff = c(1, 1.2, 0.8),
    statistic = c(3, 4, NA),
    p.value = c(0.01, 0.02, NA),
    avgAbd = 20
  )
  res <- .summary_ropeca_median_p_scaled(peptides, estimate = "diff", p.value = "p.value")
  expect_equal(res$n, 3)
  expect_equal(res$n_not_na, 2)
  # scaled p-values 0.99 and 0.98, median 0.985; the median of 2 p-values is Beta(1.5, 1.5)
  expect_equal(res$median.p.value, 0.015)
  expect_equal(res$beta.based.significance, stats::pbeta(0.015, 1.5, 1.5))
})

test_that("summary_ropeca_median_p_scaled matches a hand computation", {
  peptides <- data.frame(
    contrast = "c1",
    protein_Id = "P1",
    diff = c(1, 2, -0.5),
    statistic = c(4, 3, -1),
    p.value = c(0.01, 0.04, 0.5),
    avgAbd = c(18, 20, 22)
  )
  res <- .summary_ropeca_median_p_scaled(peptides, estimate = "diff", p.value = "p.value")
  # signed scaled p-values 1 - 0.01, 1 - 0.04 and 0.5 - 1; their median is 0.96
  expect_equal(res$median.p.scaled, 0.96)
  expect_equal(res$median.p.value, 0.04)
  # the median of 3 uniform p-values is Beta(2, 2), CDF 3 x^2 - 2 x^3
  expect_equal(res$beta.based.significance, 3 * 0.04^2 - 2 * 0.04^3)
  expect_equal(res$estimate, 1)
  expect_equal(res$mad.estimate, 1.4826)
  expect_equal(res$statistic, 3)
  expect_equal(res$avgAbd, 20)
})

test_that("summary_ropeca_median_p_scaled is conservative for uniform p-values with random signs", {
  # ropeca_nested is conservative under the null by construction; ropecaV2_nested is the calibrated variant
  withr::local_seed(11)
  n_proteins <- 4000
  peptides <- data.frame(
    contrast = "c1",
    protein_Id = rep(seq_len(n_proteins), each = 10),
    diff = stats::rnorm(n_proteins * 10),
    statistic = 0,
    p.value = stats::runif(n_proteins * 10),
    avgAbd = 1
  )
  res <- .summary_ropeca_median_p_scaled(peptides, estimate = "diff", p.value = "p.value")
  expect_lt(mean(res$beta.based.significance < 0.05), 0.01)
})

test_that("ropeca_nested equals per-peptide lm fits summarised by median and beta p-value", {
  istar <- sim_lfq_data_peptide_config(n_proteins = 15, with_missing = FALSE, seed = 4)
  lfq <- LFQData$new(istar$data, istar$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  ct <- ContrastsROPECANestedFacade$new(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"))$get_contrasts()

  # reference: one stats::lm per peptide; group A is the baseline level, so A - Ctrl is minus group_Ctrl
  data <- lfq$data_long()
  peptide_fits <- lapply(split(data, data$peptide_Id), function(d) {
    coefs <- summary(stats::lm(log2_abundance ~ group_, data = d))$coefficients
    data.frame(protein_Id = d$protein_Id[1], diff = -coefs["group_Ctrl", "Estimate"], p = coefs["group_Ctrl", 4])
  })
  peptide_fits <- do.call(rbind, peptide_fits)
  reference <- do.call(
    rbind,
    lapply(split(peptide_fits, peptide_fits$protein_Id), function(p) {
      scaled <- ifelse(p$diff > 0, 1 - p$p, p$p - 1)
      median_p <- 1 - abs(stats::median(scaled))
      shape <- min(nrow(p), 10) / 2 + 0.5
      data.frame(
        protein_Id = p$protein_Id[1],
        diff = stats::median(p$diff),
        p.value = stats::pbeta(median_p, shape, shape)
      )
    })
  )
  ct <- ct[match(reference$protein_Id, ct$protein_Id), ]
  expect_equal(ct$diff, reference$diff, tolerance = 1e-8)
  expect_equal(ct$p.value, reference$p.value, tolerance = 1e-8)
})

test_that("ContrastsTable (passive container)", {
  # Build a contrast result to feed into ContrastsTable
  istar <- sim_lfq_data_protein_config(n_proteins = 20, weight_missing = 0.4)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  fa <- ContrastsLMImputeFacade$new(lProt, "~ group_", Contr)
  ctr <- fa$get_contrasts()

  ct <- ContrastsTable$new(ctr, subject_id = fa$contrast$subject_id, model_name = "TableTest")

  x <- ct$get_contrasts()
  expect_s3_class(x, "data.frame")

  expect_null(ct$get_contrast_sides())
  expect_null(ct$get_linfct())

  tw <- ct$to_wide()
  expect_s3_class(tw, "data.frame")

  pl <- ct$get_Plotter()
  expect_s3_class(pl, "ContrastsPlotter")
})

test_that("ContrastsTable plotter uses available score columns", {
  ctr <- data.frame(
    protein_Id = paste0("P", seq_len(4)),
    contrast = "A_vs_B",
    diff = c(-2, -0.5, 0.5, 2),
    FDR = c(0.01, 0.2, 0.3, 0.02),
    modelName = "TableTest"
  )

  ct <- ContrastsTable$new(ctr, subject_id = "protein_Id", model_name = "TableTest")
  pl <- ct$get_Plotter()

  volcano_scores <- vapply(pl$volcano_spec, `[[`, character(1), "score")
  histogram_scores <- vapply(pl$histogram_spec, `[[`, character(1), "score")

  expect_equal(volcano_scores, "FDR")
  expect_equal(histogram_scores, "FDR")
  expect_named(pl$volcano(), "FDR")
  expect_s3_class(pl$volcano()$FDR, "ggplot")
})

test_that("ContrastsTable provides rank and ORA input tables", {
  ctr <- data.frame(
    protein_Id = paste0("P", seq_len(6)),
    contrast = c(rep("A_vs_B", 3), rep("C_vs_D", 3)),
    diff = c(1.5, -1.4, 0.2, 2.1, -2.2, 1.2),
    FDR = c(0.01, 0.02, 0.01, 0.2, 0.03, 0.04),
    statistic = c(5, -4, 1, 3, -6, 7),
    modelName = "TableTest"
  )

  ct <- ContrastsTable$new(ctr, subject_id = "protein_Id", model_name = "TableTest")

  rank_table <- ct$get_rank(score = "statistic")
  expect_named(rank_table, c("protein_Id", "contrast", "score"))
  expect_equal(rank_table$score, ctr$statistic)

  # Default rank score: falls back to the effect column (diff) when no p.value
  # column is present in the contrast table.
  rank_default <- ct$get_rank()
  expect_equal(rank_default$score, ctr$diff)

  ora_up <- ct$get_ora(up = TRUE, FDR_threshold = 0.05, diff_threshold = 1)
  expect_equal(ora_up$protein_Id, c("P1", "P6"))

  ora_down <- ct$get_ora(up = FALSE, FDR_threshold = 0.05, diff_threshold = 1)
  expect_equal(ora_down$protein_Id, c("P2", "P5"))
})

test_that("merge_contrasts_results labels rows when both inputs share a model name", {
  prefer_df <- data.frame(
    protein_Id = c("P1", "P2"),
    contrast = "A_vs_B",
    modelName = "lm",
    diff = c(1, 2),
    p.value = c(0.01, 0.2),
    FDR = c(0.02, 0.2),
    statistic = c(3, NA)
  )
  add_df <- transform(prefer_df, diff = c(1.1, 2.1), statistic = c(2.8, 2.5))
  prefer <- ContrastsTable$new(prefer_df, subject_id = "protein_Id", model_name = "lm")
  add <- ContrastsTable$new(add_df, subject_id = "protein_Id", model_name = "lm")

  merged <- merge_contrasts_results(prefer, add)$merged
  res <- merged$get_contrasts()
  expect_equal(as.character(res$modelName), c("lm_prefer", "lm_add"))
  expect_equal(res$diff, c(1, 2.1))
  expect_equal(merged$model_name, "lm_prefer_lm_add")
})

test_that("ContrastsFirth and ContrastsFirthFacade", {
  skip_if_not_installed("logistf")
  istar <- sim_lfq_data_protein_config(n_proteins = 20, with_missing = TRUE, weight_missing = 0.5, seed = 9)
  lfqdata <- LFQData$new(istar$data, istar$config)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  mod <- build_model_glm_protein(lfqdata, "~ group_")
  ctr <- ContrastsFirth$new(mod, Contr)
  expect_true(inherits(ctr, "ContrastsFirth"))

  res <- ctr$get_contrasts()
  expect_s3_class(res, "data.frame")
  expect_true(all(c("contrast", "diff", "statistic", "p.value", "FDR") %in% colnames(res)))

  tw <- ctr$to_wide()
  expect_s3_class(tw, "data.frame")

  pl <- ctr$get_Plotter()
  expect_s3_class(pl, "ContrastsPlotter")

  fa <- build_contrast_analysis(lfqdata, "~ group_", Contr, method = "firth")
  expect_true(inherits(fa, "ContrastsFirthFacade"))
  fa_res <- fa$get_contrasts()
  expect_false("facade" %in% colnames(fa_res))
  expect_true(all(fa_res$modelName == "firth"))
})

test_that("ContrastsFirth$get_linfct does not mutate the shared model", {
  skip_if_not_installed("logistf")
  istar <- sim_lfq_data_protein_config(n_proteins = 20, with_missing = TRUE, weight_missing = 0.5, seed = 9)
  lfqdata <- LFQData$new(istar$data, istar$config)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
  mod <- build_model_glm_protein(lfqdata, "~ group_")

  before_rows <- nrow(mod$models$model_df)
  expect_false("linfct" %in% colnames(mod$models$model_df))

  # computing contrasts must not filter rows or add columns on the shared model
  invisible(ContrastsFirth$new(mod, Contr)$get_contrasts())
  expect_equal(nrow(mod$models$model_df), before_rows)
  expect_false("linfct" %in% colnames(mod$models$model_df))

  # get_linfct returns the annotated copies, leaving the input untouched
  linfct_models <- ContrastsFirth$new(mod, Contr)$get_linfct()
  expect_true("linfct" %in% colnames(linfct_models$model_df))
  expect_false("linfct" %in% colnames(mod$models$model_df))
})

test_that("ContrastsFirth drops failed logistf fits before contrasts", {
  skip_if_not_installed("logistf")
  istar <- sim_lfq_data_protein_config(n_proteins = 20, with_missing = TRUE, weight_missing = 0.5, seed = 9)
  lfqdata <- LFQData$new(istar$data, istar$config)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  mod <- build_model_glm_protein(lfqdata, "~ group_")
  model_df <- mod$models$model_df
  failed_idx <- which(model_df$has_model_fit)[1]
  mod$models$model_df$linear_model[[failed_idx]] <- "failed fit"
  mod$models$model_df$has_model_fit[[failed_idx]] <- FALSE
  mod$models$model_df$nr_coef_not_NA[[failed_idx]] <- NA_integer_

  ctr <- ContrastsFirth$new(mod, Contr)
  res <- ctr$get_contrasts()

  expect_s3_class(res, "data.frame")
  expect_true(all(c("protein_Id", "contrast", "diff", "p.value", "FDR") %in% colnames(res)))
  expect_false(any(is.na(res$protein_Id)))
})

test_that("compute_lmer_contrast reports an empty contrast matrix instead of failing", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  # A protein whose model carries none of the requested contrast coefficients
  # leaves `linfct_matrix_contrasts()` with a 0-row matrix. lmerTest::contest()
  # indexes 1:nrow(L) and used to abort the whole run on it.
  mb <- muffle_warnings_matching(.sim_make_model_lmer("interaction"), single_peptide_warning)
  linfct <- .linfct_from_model(mb)$linfct_factors
  empty <- linfct[0, , drop = FALSE]

  expect_warning(res <- .compute_lmer_contrast(mb, empty), "No estimable contrast")
  expect_true(is.logical(res))
  expect_true(is.na(res))
})

test_that("merge_contrasts_results errors on duplicated subject-contrast rows instead of multiplying them (review B3)", {
  prefer_df <- data.frame(
    protein_Id = c("P1", "P2"),
    contrast = "A_vs_B",
    modelName = "prefer",
    diff = c(1, 2),
    p.value = c(0.01, 0.2),
    FDR = c(0.02, 0.2),
    statistic = c(3, NA)
  )
  add_df <- data.frame(
    protein_Id = c("P1", "P1", "P2"),
    contrast = "A_vs_B",
    modelName = "add",
    diff = c(1.1, 1.2, 2.1),
    p.value = c(0.02, 0.02, 0.03),
    FDR = c(0.03, 0.03, 0.04),
    statistic = c(2.8, 2.7, 2.5)
  )
  prefer <- ContrastsTable$new(prefer_df, subject_id = "protein_Id", model_name = "prefer")
  add <- ContrastsTable$new(add_df, subject_id = "protein_Id", model_name = "add")
  expect_error(merge_contrasts_results(prefer, add), "must match at most 1 row")
})

test_that("get_rank scores by signed -log10 p-value; filter_significant is one-sided for directional backends", {
  df <- data.frame(
    protein_Id = c("P1", "P2", "P3", "P4"),
    contrast = "c1",
    modelName = "m",
    diff = c(2, -2, 0.5, NA),
    statistic = c(5, -5, 1, NA),
    p.value = c(0.001, 0.01, 0.5, NA),
    FDR = c(0.004, 0.02, 0.5, NA),
    avgAbd = 1
  )
  ct <- ContrastsTable$new(df, subject_id = "protein_Id")
  rank <- ct$get_rank()
  expect_identical(rank$protein_Id, c("P1", "P2", "P3"))
  expect_equal(rank$score, c(3, -2, -log10(0.5)))
  expect_identical(ct$filter_significant(FDR_threshold = 0.05, diff_threshold = 1)$protein_Id, c("P1", "P2"))

  directional <- ContrastsTable$new(
    df,
    subject_id = "protein_Id",
    config = ContrastConfiguration$new(subject_id = "protein_Id", significance_directional = TRUE)
  )
  expect_identical(directional$filter_significant(FDR_threshold = 0.05, diff_threshold = 1)$protein_Id, "P1")
})
