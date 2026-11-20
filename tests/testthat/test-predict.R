# predict() on fitted analyses: new combinations of subject and group (srr RE4.16).

predict_fixture <- function() {
  sim <- sim_lfq_data_protein_config(Nprot = 6, with_missing = FALSE, seed = 1)
  LFQData$new(sim$data, sim$config)
}

join_fitted <- function(prediction, fitted) {
  by <- intersect(intersect(names(prediction), names(fitted)), c("protein_Id", "peptide_Id", "sampleName", "group_"))
  dplyr::inner_join(prediction, fitted, by = by)
}

test_that("predict() reproduces fitted values on observed cells", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- predict_fixture()
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  for (method in c("lm", "rlm", "rfit", "limma", "deqms", "firth", "lm_impute")) {
    if (method == "rfit") {
      skip_if_not_installed("Rfit")
    }
    fit <- suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrasts, method))
    prediction <- predict(fit)
    expect_identical(nrow(prediction), 6L * 12L, label = method)
    expect_named(prediction, c("protein_Id", "sample", "sampleName", "group_", "prediction", "imputed"))
    joined <- join_fitted(prediction, fitted(fit))
    expect_equal(joined$prediction, joined$fitted, tolerance = 1e-8, label = method)
  }
})

test_that("nested methods predict per subject of their model", {
  withr::local_options(prolfqua.verbose = FALSE)
  sim <- sim_lfq_data_peptide_config(Nprot = 5, with_missing = FALSE, seed = 1)
  peptide <- suppressMessages(LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq)
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  ropeca <- build_contrast_analysis(peptide, "~ group_", contrasts, "ropeca_nested")
  joined <- join_fitted(predict(ropeca), fitted(ropeca))
  expect_true("peptide_Id" %in% names(predict(ropeca)))
  expect_equal(joined$prediction, joined$fitted, tolerance = 1e-8)
  binomial <- suppressWarnings(build_contrast_analysis(peptide, "~ group_", contrasts, "binomial_nested"))
  probabilities <- predict(binomial)$prediction
  expect_true(all(probabilities >= 0 & probabilities <= 1))
  hurdle <- suppressWarnings(build_contrast_analysis(peptide, "~ group_", contrasts, "hurdle_nested"))
  expect_error(predict(hurdle), class = "prolfqua_error_not_available")
  expect_s3_class(predict(hurdle$contrast$intensity), "data.frame")
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  lmer <- suppressMessages(suppressWarnings(
    build_contrast_analysis(peptide, "~ group_ + (1 | peptide_Id)", contrasts, "lmer_nested")
  ))
  prediction <- predict(lmer)
  fitted_subjects <- lmer$model$model_df$protein_Id[lmer$model$model_df$has_model_fit]
  expect_false(anyNA(prediction$prediction[prediction$protein_Id %in% fitted_subjects]))
  expect_true(all(is.na(prediction$prediction[!prediction$protein_Id %in% fitted_subjects])))
})

test_that("a group absent from one protein is predicted from the LOD refit, NA without it", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- predict_fixture()
  absent <- unique(lfq$data_long()$protein_Id)[1]
  data <- lfq$data_long()
  data$abundance[data$protein_Id == absent & data$group_ == "A"] <- NA
  lfq$set_data(data)
  contrasts <- c(AvsC = "group_A - group_Ctrl")
  newdata <- data.frame(protein_Id = absent, group_ = c("A", "Ctrl"))
  plain <- predict(suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrasts, "lm")), newdata)
  expect_true(is.na(plain$prediction[1]))
  expect_false(is.na(plain$prediction[2]))
  refit <- predict(suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrasts, "lm_impute")), newdata)
  expect_false(anyNA(refit$prediction))
  expect_true(all(refit$imputed))
  limma <- predict(suppressWarnings(build_contrast_analysis(lfq, "~ group_", contrasts, "limma_impute")), newdata)
  expect_false(anyNA(limma$prediction))
})

test_that("unknown groups and subjects give NA; newdata must name the subject", {
  withr::local_options(prolfqua.verbose = FALSE)
  lfq <- predict_fixture()
  fit <- build_contrast_analysis(lfq, "~ group_", c(AvsC = "group_A - group_Ctrl"), "lm")
  subject <- unique(lfq$data_long()$protein_Id)[1]
  newdata <- data.frame(protein_Id = c(subject, subject, "unknown"), group_ = c("A", "Z", "A"))
  expect_identical(is.na(predict(fit, newdata)$prediction), c(FALSE, TRUE, TRUE))
  expect_error(predict(fit, data.frame(group_ = "A")), class = "prolfqua_error")
  expect_error(predict(fit, list(protein_Id = subject)), class = "prolfqua_error_bad_argument")
})
