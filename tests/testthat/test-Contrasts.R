# Test Contrasts R6 classes

test_that("Contrasts (Wald test)", {
  istar <- sim_lfq_data_peptide_config(Nprot = 20)
  config <- istar$config

  modelFunction <-
    strategy_lmer("abundance ~ group_ + (1 | peptide_Id) + (1 | sample)")

  mod <- build_model(
    istar$data,
    modelFunction,
    subject_Id = config$hierarchy_keys_depth()
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
    subject_Id = config$hierarchy_keys_depth()
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

  ct <- ContrastsTable$new(ctr, subject_Id = csi$subject_Id, model_name = "TableTest")

  x <- ct$get_contrasts()
  expect_s3_class(x, "data.frame")

  expect_null(ct$get_contrast_sides())
  expect_null(ct$get_linfct())

  tw <- ct$to_wide()
  expect_s3_class(tw, "data.frame")

  pl <- ct$get_Plotter()
  expect_s3_class(pl, "ContrastsPlotter")
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
