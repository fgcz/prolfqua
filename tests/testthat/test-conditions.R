#' Statistical condition checks
#'
#' @srrstatsTODO {G5.2} Condition tests exercise validation and warning paths; test-power_fdr.R covers all
#'   planning-input errors and test-warnings.R covers raw pilot input; limit: informative-message paths lack tests.
#' @srrstatsTODO {G5.2b} Coverage-guided condition tests trigger new checks; limit: informative-message paths
#'   remain.
#' @noRd

# Typed error conditions raised at prolfqua's public boundaries.

test_that("setup_analysis aborts (typed) when the file_name column is missing", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  cfg <- istar$config$clone(deep = TRUE)
  cfg$file_name <- "column_that_does_not_exist"

  expect_error(
    prolfqua::setup_analysis(istar$data, cfg),
    class = "prolfqua_error_missing_column"
  )
})

test_that("setup_analysis aborts (typed) when no factors are configured", {
  raw <- data.frame(
    sample = c("s1", "s2"),
    proteinID = c("P1", "P1"),
    idtype2 = c("x", "x"),
    peptideID = c("pep1", "pep1"),
    abundance = c(10, 12),
    stringsAsFactors = FALSE
  )
  cfg <- prolfqua::AnalysisConfiguration$new()
  cfg$file_name <- "sample"
  cfg$hierarchy[["protein_Id"]] <- c("proteinID", "idtype2")
  cfg$hierarchy[["peptide_Id"]] <- "peptideID"
  cfg$set_response("abundance")

  expect_error(
    prolfqua::setup_analysis(raw, cfg),
    class = "prolfqua_error_invalid_configuration"
  )
})

test_that("setup_analysis aborts (typed) on duplicate hierarchy-key / sample combinations", {
  raw <- data.frame(
    sample = c("s1", "s1"),
    group = c("A", "A"),
    proteinID = c("P1", "P1"),
    idtype2 = c("x", "x"),
    peptideID = c("pep1", "pep1"),
    abundance = c(10, 12),
    isotopeLabel = "light",
    qValue = 0,
    stringsAsFactors = FALSE
  )
  cfg <- prolfqua::AnalysisConfiguration$new()
  cfg$file_name <- "sample"
  cfg$factors["group_"] <- "group"
  cfg$hierarchy[["protein_Id"]] <- c("proteinID", "idtype2")
  cfg$hierarchy[["peptide_Id"]] <- "peptideID"
  cfg$set_response("abundance")

  # raw has no nr_children column, which setup_analysis() warns about before it checks the keys
  expect_error(
    muffle_warnings_matching(prolfqua::setup_analysis(raw, cfg), "no nr_children column"),
    class = "prolfqua_error_duplicate_keys"
  )
})

test_that("LFQData$set_data rejects a non-data-frame (typed)", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_error(
    lfq$set_data(list(a = 1)),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("LFQData$set_data enforces the configured schema (typed)", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_error(
    lfq$set_data(data.frame(unrelated = 1)),
    class = "prolfqua_error_missing_column"
  )
})

test_that("build_contrast_analysis aborts (typed) on an unknown method", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  lfq$rename_response("transformedIntensity")
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  expect_error(
    prolfqua::build_contrast_analysis(lfq, "~ group_", contrasts, method = "definitely_not_real"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("LFQData$get_Aggregator aborts (typed) on an unknown method", {
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_error(
    lfq$get_Aggregator("not_a_method"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("aggregating single-level hierarchy aborts (typed)", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_error(
    lfq$get_Aggregator("medpolish"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("center_to_reference aborts (typed) on mismatched transformation state", {
  istar <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  ref <- lfq$get_copy()
  ref <- ref$get_Transformer()$log2()$lfq

  expect_error(
    lfq$get_Transformer()$center_to_reference(ref),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("informative messages carry class prolfqua_message", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)

  expect_message(
    lfq$get_Transformer()$log2(),
    "Column added",
    class = "prolfqua_message"
  )
})

test_that("options(prolfqua.verbose = FALSE) silences informative messages", {
  istar <- prolfqua::sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- prolfqua::LFQData$new(istar$data, istar$config)
  withr::local_options(prolfqua.verbose = FALSE)

  expect_no_message(lfq$get_Transformer()$log2())
})

test_that("errors carry class prolfqua_error", {
  expect_error(
    prolfqua::register_facade(name = 1, class = "X", needs = "aggregated"),
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    prolfqua::contrasts_linfct(data.frame(protein_Id = c("p1", "p2")), linfct = list(1)),
    class = "prolfqua_error"
  )
})

test_that("plot_pca rejects a PC argument of the wrong length", {
  withr::local_seed(1)
  x <- matrix(rnorm(20), nrow = 5, dimnames = list(NULL, paste0("s", 1:4)))
  annot <- data.frame(sample = paste0("s", 1:4), group = c("a", "a", "b", "b"))
  expect_error(
    prolfqua::plot_pca(x, annot, sample_name = "sample", factor_keys = "group", PC = 1),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("LFQData$data_wide rejects an unknown response value", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  expect_error(LFQData$new(x$data, x$config)$data_wide(value = "absent"), class = "prolfqua_error_bad_argument")
})

test_that("aggregators reject input at the deepest hierarchy level", {
  x <- sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
  x$config$hierarchy_depth <- 2
  expect_error(LFQData$new(x$data, x$config)$get_Aggregator(), class = "prolfqua_error_bad_argument")
})

test_that("aggregator plotting requires an aggregation result", {
  x <- sim_lfq_data_peptide_config(n_proteins = 5, seed = 42)
  lfq <- LFQData$new(x$data, x$config)$get_Transformer()$log2()$lfq
  expect_error(lfq$get_Aggregator()$plot(), class = "prolfqua_error_bad_argument")
})

test_that("setup requires a file-name configuration", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  x$config$file_name <- NULL
  expect_error(setup_analysis(x$data, x$config), class = "prolfqua_error_invalid_configuration")
})

test_that("wide values must be a named list of matrices", {
  expect_error(.check_wide_values(list(1)), class = "prolfqua_error_bad_argument")
})

test_that("wide matrices must have matching shapes and samples", {
  expect_error(
    .check_wide_values(list(a = matrix(1, 2, 2), b = matrix(1, 2, 3))),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("wide rowdata must carry configured hierarchy columns", {
  cfg <- AnalysisConfiguration$new()
  cfg$hierarchy <- list(protein = "protein")
  expect_error(
    .check_wide_input(list(a = matrix(1)), data.frame(x = 1), data.frame(), cfg),
    class = "prolfqua_error_missing_column"
  )
})

test_that("wide annotation must carry the sample-name column", {
  cfg <- AnalysisConfiguration$new()
  cfg$hierarchy <- list(protein = "protein")
  expect_error(
    .check_wide_input(list(a = matrix(1)), data.frame(protein = "p1"), data.frame(), cfg),
    class = "prolfqua_error_missing_column"
  )
})

test_that("wide rowdata must match the matrix row count", {
  cfg <- AnalysisConfiguration$new()
  cfg$hierarchy <- list(protein = "protein")
  expect_error(
    .check_wide_input(list(a = matrix(1, 2, 1)), data.frame(protein = "p1"), data.frame(sampleName = "s1"), cfg),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("wide annotation must match matrix column names uniquely", {
  cfg <- AnalysisConfiguration$new()
  cfg$hierarchy <- list(protein = "protein")
  expect_error(
    .check_wide_input(
      list(a = matrix(1, dimnames = list(NULL, "s1"))),
      data.frame(protein = "p1"),
      data.frame(sampleName = "s2"),
      cfg
    ),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("contrast configuration rejects unknown fields", {
  expect_error(list_to_contrast_configuration(list(unknown_field = 1)), class = "prolfqua_error")
})

test_that("abstract model methods fail with classed errors", {
  MinimalModel <- R6::R6Class("MinimalModel", inherit = ModelInterface)
  model <- MinimalModel$new()
  expect_error(model$get_coefficients(), class = "prolfqua_error")
  expect_error(model$get_anova(), class = "prolfqua_error")
})

test_that("abstract contrast methods fail with classed errors", {
  MinimalContrast <- R6::R6Class("MinimalContrast", inherit = ContrastsInterface)
  contrast <- MinimalContrast$new()
  expect_error(contrast$get_contrast_sides(), class = "prolfqua_error")
  expect_error(contrast$get_contrasts(), class = "prolfqua_error")
  expect_error(contrast$get_missing(), class = "prolfqua_error")
})

test_that("abstract strategy fitting fails with a classed error", {
  MinimalStrategy <- R6::R6Class(
    "MinimalStrategy",
    inherit = StrategyBase,
    public = list(trigger_fit = function(x) private$fit(x))
  )
  expect_error(MinimalStrategy$new("y ~ x")$trigger_fit(data.frame()), class = "prolfqua_error")
})

test_that("Firth model builder checks class before hierarchy access", {
  expect_error(build_model_logistf(list(), "y ~ x"), class = "prolfqua_error_bad_argument")
})

test_that("detection preprocessing requires LFQData", {
  expect_error(.prepare_detection_lfqdata(list()), class = "prolfqua_error_bad_argument")
})

test_that("facade registration rejects invalid class strings", {
  expect_error(register_facade("x", class = "", needs = "aggregated"), class = "prolfqua_error_bad_argument")
})

test_that("facade registration rejects non-string input shapes", {
  expect_error(register_facade("x", class = "X", needs = NA_character_), class = "prolfqua_error_bad_argument")
})

test_that("facade registration rejects unknown input shapes", {
  expect_error(register_facade("x", class = "X", needs = "same"), class = "prolfqua_error_bad_argument")
})

test_that("facade registration rejects non-string package names", {
  expect_error(
    register_facade("x", class = "X", needs = "aggregated", package = 1),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("limma borrowing fails when no variance donor is usable", {
  fit <- list(coefficients = matrix(NA_real_, 2, 2), sigma = c(0, 0), df.residual = c(0, 0))
  expect_error(compute_borrowed_variance_limma(fit), class = "prolfqua_error")
})

test_that("model borrowing fails when no fitted donor is usable", {
  expect_error(
    compute_borrowed_variance(data.frame(has_model_fit = TRUE, isSingular = TRUE, nr_coef_not_NA = 1, df.residual = 3)),
    class = "prolfqua_error"
  )
})

test_that("contrast calculation rejects an all-zero coefficient vector", {
  model <- stats::lm(Sepal.Length ~ Species, data = iris)
  zero <- matrix(0, nrow = 1, ncol = length(stats::coef(model)), dimnames = list("zero", names(stats::coef(model))))
  expect_error(
    compute_contrast(model, zero),
    regexp = "`linfct` must have a non-zero coefficient in every row",
    class = "prolfqua_error"
  )
})

test_that("imputation rejects different subject identifiers", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- LFQData$new(x$data, x$config)
  model <- build_model(lfq, strategy_lm("abundance ~ group_"))
  model$subject_id <- "other"
  expect_error(impute_from_model(model, lfq), class = "prolfqua_error")
})

test_that("imputation rejects non-lm fits", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- LFQData$new(x$data, x$config)
  model <- build_model(lfq, strategy_lm("abundance ~ group_"))
  model$model_df$linear_model[[1]] <- list(unsupported = TRUE)
  expect_error(impute_from_model(model, lfq), class = "prolfqua_error")
})

test_that("imputation rejects mismatched responses", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- LFQData$new(x$data, x$config)
  model <- build_model(lfq, strategy_lm("abundance ~ group_"))
  lfq$rename_response("other")
  expect_error(impute_from_model(model, lfq), class = "prolfqua_error")
})

test_that("imputation requires one fitted-model row for every input subject", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  lfq <- LFQData$new(x$data, x$config)
  model <- build_model(lfq, strategy_lm("abundance ~ group_"))
  model$model_df <- model$model_df[-1, ]
  expect_error(impute_from_model(model, lfq), class = "prolfqua_error")
})

test_that("unknown progress reporters carry the package error class", {
  expect_error(.make_progress(3, reporter = "unknown"), class = "prolfqua_error")
})

test_that("limma rejects unsupported weights specifications", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  expect_error(
    build_model_limma(LFQData$new(x$data, x$config), strategy_limma("abundance ~ group_", weights = list(1))),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("Hurdle requires a contrast interface for intensity", {
  expect_error(ContrastsHurdleFacade$new(list(), ContrastsInterface$new()), class = "prolfqua_error_bad_argument")
})

test_that("Hurdle requires a contrast interface for detection", {
  expect_error(ContrastsHurdleFacade$new(ContrastsInterface$new(), list()), class = "prolfqua_error_bad_argument")
})

test_that("facades check LFQData even when called without the factory", {
  expect_error(ContrastsLMFacade$new(list(), "~ group_", c(a = "A - B")), class = "prolfqua_error_bad_argument")
})

test_that("detection count summaries reject non-LFQData input", {
  expect_error(.summarize_detection_counts(list()), class = "prolfqua_error_bad_argument")
})

test_that("detection count summaries reject negative detection counts", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  x$data$bin_resp <- -1
  x$config$bin_resp <- "bin_resp"
  expect_error(.summarize_detection_counts(LFQData$new(x$data, x$config)), class = "prolfqua_error")
})

test_that("plot label widths reject truncation below the required minimum", {
  expect_error(.plot_label_widths("label", max_chars = 1, min_chars = 3), class = "prolfqua_error_bad_argument")
})

test_that("NA heatmap annotation cannot name unknown samples", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  expect_error(
    plot_na_heatmap(matrix(1, dimnames = list(NULL, "s1")), data.frame(sample = "s2"), "group", "sample"),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("PCA rejects a requested component beyond the available rank", {
  x <- matrix(seq_len(8), nrow = 2, dimnames = list(NULL, paste0("s", 1:4)))
  annot <- data.frame(sample = paste0("s", 1:4), group = "A")
  expect_error(plot_pca(x, annot, "sample", "group", PC = c(1, 3)), class = "prolfqua_error")
})

test_that("matrix conversion rejects other object types", {
  expect_error(matrix_to_tibble(list()), class = "prolfqua_error_bad_argument")
})

test_that("missing contrast columns raise a classed error", {
  expect_error(
    .stop_if_missing_columns(data.frame(), "p.value", "filter contrasts"),
    class = "prolfqua_error_missing_column"
  )
})

test_that("original ROPECA rejects aggregated contrasts", {
  expect_error(
    ContrastsROPECA$new(ContrastsTable$new(data.frame(protein_Id = "p1"))),
    class = "prolfqua_error_bad_argument"
  )
})

test_that("nested mixed models require a child hierarchy for automatic random effects", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  expect_error(.augment_lmer_nested_modelstr("~ group_", LFQData$new(x$data, x$config)), class = "prolfqua_error")
})

test_that("limpa requires a configured child count column", {
  skip_if_not_installed("limpa")
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  x$config$nr_children <- NULL
  expect_error(
    build_model_limpa(LFQData$new(x$data, x$config), strategy_limpa("abundance ~ group_")),
    class = "prolfqua_error_invalid_configuration"
  )
})

test_that("error helpers name the calling function in the header, not themselves (review B2)", {
  user_fn <- function(x) prolfqua:::abort_bad_argument("x", "be numeric", not = "a string")
  cnd <- rlang::catch_cnd(user_fn(1))
  expect_s3_class(cnd, "prolfqua_error_bad_argument")
  expect_identical(rlang::call_name(cnd$call), "user_fn")
  expect_match(conditionMessage(cnd), "`x` must be numeric.", fixed = TRUE)
  expect_match(conditionMessage(cnd), "You supplied a string.", fixed = TRUE)

  user_cols <- function(df) prolfqua:::abort_missing_columns(c("a", "b"), "df")
  cnd <- rlang::catch_cnd(user_cols(1))
  expect_s3_class(cnd, "prolfqua_error_missing_column")
  expect_identical(rlang::call_name(cnd$call), "user_cols")
  expect_match(conditionMessage(cnd), "Can't find columns `a`, `b` in `df`.", fixed = TRUE)

  user_config <- function(config) prolfqua:::abort_invalid_config("`config` must name a sample column.")
  cnd <- rlang::catch_cnd(user_config(1))
  expect_s3_class(cnd, "prolfqua_error_invalid_configuration")
  expect_identical(rlang::call_name(cnd$call), "user_config")

  # internal checkers forward their caller, also from inside a tryCatch() handler
  user_formula <- function(modelstr) prolfqua:::.as_model_formula(modelstr)
  expect_identical(rlang::call_name(rlang::catch_cnd(user_formula(1))$call), "user_formula")
  user_choice <- function(method) prolfqua:::.arg_match_prolfqua(method, c("lm", "limma"), "method")
  cnd <- rlang::catch_cnd(user_choice("lmm"))
  expect_identical(rlang::call_name(cnd$call), "user_choice")
  expect_match(conditionMessage(cnd), "`method` must be one of \"lm\", \"limma\".", fixed = TRUE)

  cnd <- rlang::catch_cnd(sim_lfq_data_peptide_config(n_proteins = -1))
  expect_identical(rlang::call_name(cnd$call), "sim_lfq_data_peptide_config")
})
