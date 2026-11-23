#' Statistical condition checks
#'
#' @srrstats {G5.2a} Message texts are unique: informative messages and warnings name their operation. Limit: an
#'   argument check reused across functions (e.g. "`lfqdata` must be an LFQData object") keeps one wording; the error
#'   header names the function that was called.
#' @noRd

test_that("log transforms set zeros to NA and return which", {
  for (fun in list(log, log2, log10)) {
    res <- .transform_work_intensity(data.frame(x = c(0, 1)), "x", fun, "f")
    expect_identical(res$data$f_x, c(NA, 0))
    expect_identical(res$zeros, c(TRUE, FALSE))
  }
})

test_that("log transforms give one classed warning when negative values become NaN", {
  expect_warning(
    res <- .transform_work_intensity(data.frame(x = c(-1, 1, 0)), "x", log2),
    "1 negative values.*become NaN",
    class = "prolfqua_negative_intensity"
  )
  expect_identical(res$data$log2_x, c(NaN, 0, NA))
  expect_identical(res$data$x, c(-1, 1, NA))
})

test_that("repeated log transformation warns", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  trans <- LFQData$new(x$data, x$config)$get_Transformer()$log2()
  expect_warning(trans$log2(), "already transformed")
})

test_that("zero MAD warns and skips division by zero", {
  expect_warning(result <- .robust_scale(matrix(1, 3, 3)), "MAD values are zero")
  expect_true(all(is.finite(result)))
})

test_that("centering raw data warns about the required log scale", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- LFQData$new(x$data, x$config)
  expect_warning(lfq$get_Transformer()$center_to_reference(lfq), "should be log2-transformed")
})

test_that("median polish warns on raw intensities", {
  x <- sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
  expect_warning(LFQData$new(x$data, x$config)$get_Aggregator("medpolish"), "medpolish works best")
})

test_that("robust aggregation warns on raw intensities", {
  x <- sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
  expect_warning(LFQData$new(x$data, x$config)$get_Aggregator("rlm"), "Robust regression works best")
})

test_that("top-N aggregation warns on transformed intensities", {
  x <- sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
  lfq <- LFQData$new(x$data, x$config)$get_Transformer()$log2()$lfq
  expect_warning(lfq$get_Aggregator("topN"), "top N works with raw data")
})

test_that("empty heatmaps warn with the plotting function name", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  expect_warning(
    .plot_heatmap(matrix(numeric(), 0, 2), data.frame(), "group", "sample"),
    "plot_heatmap: the dataset has no rows"
  )
})

test_that("single-row rasters warn with the plotting function name", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  expect_warning(
    .plot_raster(matrix(1, 1, 2), data.frame(), "group", "sample"),
    "plot_raster: the dataset needs at least two rows"
  )
})

test_that("Firth ANOVA reports its unimplemented method", {
  model <- ModelFirth$new(list(model_df = data.frame()))
  expect_warning(model$get_anova(), "`get_anova\\(\\)` is not implemented")
})

test_that("Firth ANOVA histogram reports its unimplemented method", {
  model <- ModelFirth$new(list(model_df = data.frame()))
  expect_warning(model$anova_histogram(), "`anova_histogram\\(\\)` is not implemented")
})

test_that("failed contrasts report how many proteins were dropped", {
  models <- tibble::tibble(protein_Id = "p1", linear_model = list(NULL), isSingular = TRUE, sigma = 1, df.residual = 3)
  expect_warning(
    .contrasts_linfct(models, matrix(1), contrastfun = function(...) NA),
    "contrasts_linfct: dropped 1 of 1 proteins"
  )
})

test_that("mixed-model contrast dimension mismatches warn", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("lme4")
  model <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  expect_warning(.compute_lmer_contrast(model, matrix(1, 1, 3)), "rank-deficient")
})

test_that("empty mixed-model contrasts warn", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("lme4")
  model <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  expect_warning(.compute_lmer_contrast(model, matrix(numeric(), 0, 2)), "No estimable contrast")
})

test_that("accession parsing warns when no accessions match", {
  expect_warning(get_uniprot_id_from_fasta_header(data.frame(protein_Id = "not-an-accession")), "No SwissProt")
})

test_that("heatmaps warn when every feature has zero variation", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  x <- matrix(1, 2, 3, dimnames = list(NULL, paste0("s", 1:3)))
  expect_warning(
    .plot_heatmap(x, data.frame(sample = paste0("s", 1:3), group = "A"), "group", "sample"),
    "No features have enough observed variation"
  )
})

test_that("MA plots warn when average abundance is unavailable", {
  table <- data.frame(protein_Id = "p1", contrast = "a", diff = 1, p.value = 0.1, FDR = 0.1, modelName = "lm")
  expect_warning(ContrastsPlotter$new(table, "protein_Id")$ma_plot(), "Can't draw an MA plot")
})

test_that("limma ANOVA warns for an intercept-only model", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  model <- build_model_limma(LFQData$new(x$data, x$config), strategy_limma("abundance ~ 1"))
  expect_warning(model$get_anova(), "No non-intercept coefficients for ANOVA")
})

test_that("Firth fitting returns the failure when coefficient iteration cannot converge", {
  dat <- data.frame(
    group_ = rep(c("A", "B"), each = 4),
    detected = c(1, 2, 1, 3, 4, 5, 3, 5),
    undetected = c(4, 3, 4, 2, 1, 0, 2, 0)
  )
  original_coef <- stats::coef
  calls <- 0L
  local_mocked_bindings(
    coef = function(object, ...) {
      calls <<- calls + 1L
      original_coef(object, ...) + calls
    },
    .package = "stats"
  )
  expect_match(strategy_binomial_firth("~ group_")$model_fun(dat), "Firth iteration did not converge")
  expect_gte(calls, 100L)
})
