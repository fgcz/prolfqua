# Test Contrasts R6 classes

test_that("Contrasts (Wald test)", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20)
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

test_that("ContrastsModerated", {
  istar <- sim_lfq_data_protein_config(Nprot = 20)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  modelFunction <- strategy_lm("transformedIntensity ~ group_")
  mod <- build_model(lProt, modelFunction)

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  contrast <- Contrasts$new(mod, Contr)
  moderated <- ContrastsModerated$new(contrast)

  x <- moderated$get_contrasts()
  expect_s3_class(x, "data.frame")
  expect_true(all(c("diff", "p.value", "FDR") %in% colnames(x)))

  tw <- moderated$to_wide()
  expect_s3_class(tw, "data.frame")

  pl <- moderated$get_Plotter()
  expect_s3_class(pl, "ContrastsPlotter")
})

test_that("ContrastsROPECA", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20)
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

test_that("ContrastsMissing", {
  istar <- sim_lfq_data_protein_config(Nprot = 20, weight_missing = 0.4)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  csi <- ContrastsMissing$new(lProt, contrasts = Contr)

  res <- csi$get_contrasts()
  expect_s3_class(res, "data.frame")
  expect_true(all(c("diff", "p.value", "FDR") %in% colnames(res)))
  expect_equal(sum(is.na(res$p.value)), 0)

  cs <- csi$get_contrast_sides()
  expect_s3_class(cs, "data.frame")

  pl <- csi$get_Plotter()
  expect_s3_class(pl, "ContrastsPlotter")
})

test_that("ContrastsTable (passive container)", {
  # Build a contrast result to feed into ContrastsTable
  istar <- sim_lfq_data_protein_config(Nprot = 20, weight_missing = 0.4)
  lProt <- LFQData$new(istar$data, istar$config)
  lProt$rename_response("transformedIntensity")

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")
  csi <- ContrastsMissing$new(lProt, contrasts = Contr)
  ctr <- csi$get_contrasts()

  ct <- ContrastsTable$new(ctr, subject_id = csi$subject_id, model_name = "TableTest")

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

test_that("ContrastsFirth and ContrastsFirthFacade", {
  istar <- sim_lfq_data_protein_config(Nprot = 20, with_missing = TRUE, weight_missing = 0.5, seed = 9)
  lfqdata <- LFQData$new(istar$data, istar$config)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  mod <- build_model_glm_protein(lfqdata, "~ group_")
  ctr <- ContrastsFirth$new(mod, Contr)

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
  expect_true("facade" %in% colnames(fa_res))
  expect_true(all(fa_res$facade == "firth"))
})

test_that("ContrastsFirth drops failed logistf fits before contrasts", {
  istar <- sim_lfq_data_protein_config(Nprot = 20, with_missing = TRUE, weight_missing = 0.5, seed = 9)
  lfqdata <- LFQData$new(istar$data, istar$config)
  Contr <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  mod <- build_model_glm_protein(lfqdata, "~ group_")
  model_df <- mod$models$models1$model_df
  failed_idx <- which(model_df$has_model_fit)[1]
  mod$models$models1$model_df$linear_model[[failed_idx]] <- "failed fit"
  mod$models$models1$model_df$has_model_fit[[failed_idx]] <- FALSE
  mod$models$models1$model_df$nr_coef_not_NA[[failed_idx]] <- NA_integer_

  ctr <- ContrastsFirth$new(mod, Contr)
  res <- ctr$get_contrasts()

  expect_s3_class(res, "data.frame")
  expect_true(all(c("protein_Id", "contrast", "diff", "p.value", "FDR") %in% colnames(res)))
  expect_false(any(is.na(res$protein_Id)))
})
