# Input types, lengths and conversions at the input boundary (srr G2.0-G2.12)

types_input <- function(factor_column = "group") {
  raw <- withr::with_seed(1, sim_lfq_data(n_proteins = 6, n_per_group = 3))
  raw$nr_children <- 1
  config <- AnalysisConfiguration$new()
  config$file_name <- "sample"
  config$factors[paste0(factor_column, "_")] <- factor_column
  config$hierarchy[["protein_Id"]] <- "proteinID"
  config$set_response("abundance")
  list(data = as.data.frame(raw), config = config)
}

test_that("single-valued arguments reject vectors, NA, wrong types and non-whole counts", {
  withr::local_options(prolfqua.verbose = FALSE)
  x <- types_input()
  bad <- function(expr, arg) {
    expect_error(expr, paste0("`", arg, "`"), class = "prolfqua_error_bad_argument")
  }
  bad(LFQData$new(x$data, x$config, prefix = 1:3), "prefix")
  bad(LFQData$new(x$data, x$config, setup = NA), "setup")
  bad(setup_analysis(x$data, x$config, debug = "no"), "debug")
  bad(sim_lfq_data(n_per_group = 2.5), "n_per_group")
  bad(sim_lfq_data(n_proteins = c(5, 6)), "n_proteins")
  bad(sim_lfq_data(probability_of_success = 2), "probability_of_success")
  bad(sim_lfq_data(mean_prot = -1), "mean_prot")
  bad(sim_lfq_data(peptide = "yes"), "peptide")
  bad(sim_lfq_data_peptide_config(seed = c(1, 2)), "seed")
  bad(sim_lfq_data_peptide_config(n_per_group = 0), "n_per_group")
  bad(sim_lfq_data_protein_config(with_missing = c(TRUE, FALSE)), "with_missing")
  bad(sim_lfq_data_protein_config(paired = NA), "paired")
  bad(sim_lfq_data_2factor_config(two_factor = NA), "two_factor")
  bad(which_missing(c(1.5, 2.5, 3.5), weight_missing = c(0.1, 0.2)), "weight_missing")
  bad(strategy_limma("~ group_", trend = NA), "trend")
  bad(strategy_limpa("~ group_", robust = "no"), "robust")
  bad(strategy_lm("~ group_", model_name = 1), "model_name")
  # counts are converted to integer
  expect_identical(.assert_count(3), 3L)
  expect_identical(.assert_count(0, min = 0L), 0L)
})

test_that("null_split_check takes whole split and feature counts and one seed", {
  withr::local_options(prolfqua.verbose = FALSE)
  sim <- sim_lfq_data_protein_config(n_proteins = 10, seed = 1)
  lfq <- LFQData$new(sim$data, sim$config)
  # 2.5 splits used to give 98 splits: the loop stopped only when the count matched exactly
  expect_error(null_split_check(lfq, n_splits = 2.5), "`n_splits`", class = "prolfqua_error_bad_argument")
  expect_error(null_split_check(lfq, n_features = -1), "`n_features`", class = "prolfqua_error_bad_argument")
  expect_error(null_split_check(lfq, seed = c(1, 2)), "`seed`", class = "prolfqua_error_bad_argument")
})

test_that("sample size planning checks delta, power, significance level and quantiles", {
  withr::local_options(prolfqua.verbose = FALSE)
  sim <- sim_lfq_data_protein_config(n_proteins = 10, seed = 1)
  stats <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq$get_Stats()
  bad <- function(expr, arg) {
    expect_error(expr, paste0("`", arg, "`"), class = "prolfqua_error_bad_argument")
  }
  bad(stats$power_t_test(sig.level = c(0.01, 0.05)), "sig.level")
  bad(stats$power_t_test(power = 1), "power")
  bad(stats$power_t_test(delta = -1), "delta")
  bad(stats$power_t_test_quantiles(probs = 2), "probs")
  expect_gt(nrow(stats$power_t_test(delta = 1)), 0)
})

test_that("a factor column keeps its level order, which sets the reference level", {
  withr::local_options(prolfqua.verbose = FALSE)
  x <- types_input()
  contrasts <- c(AvsCtrl = "group_A - group_Ctrl", BvsCtrl = "group_B - group_Ctrl")
  as_character <- build_contrast_analysis(LFQData$new(x$data, x$config, setup = TRUE), "~ group_", contrasts, "lm")

  x$data$group <- factor(x$data$group, levels = c("Ctrl", "B", "A"))
  lfq <- LFQData$new(x$data, x$config, setup = TRUE)
  expect_identical(levels(lfq$data_long()$group_), c("Ctrl", "B", "A"))
  as_factor <- build_contrast_analysis(lfq, "~ group_", contrasts, "lm")
  fit <- as_factor$model$model_df$linear_model[[1]]
  expect_identical(names(stats::coef(fit)), c("(Intercept)", "group_B", "group_A"))
  # the reference level changes the coefficients, not the contrasts
  expect_equal(as_factor$get_contrasts()$diff, as_character$get_contrasts()$diff, tolerance = 1e-10)

  withr::local_options(prolfqua.verbose = TRUE)
  x$data$group <- factor(x$data$group, levels = c("Ctrl", "B", "A", "unused"), ordered = TRUE)
  expect_message(
    lfq <- LFQData$new(x$data, x$config, setup = TRUE),
    "kept as an unordered factor with levels Ctrl, B, A",
    class = "prolfqua_message"
  )
  expect_false(is.ordered(lfq$data_long()$group_))
  expect_identical(levels(lfq$data_long()$group_), c("Ctrl", "B", "A"))
})

test_that("an integer response is stored as double", {
  withr::local_options(prolfqua.verbose = FALSE)
  x <- types_input()
  x$data$abundance <- as.integer(round(x$data$abundance))
  expect_type(setup_analysis(x$data, x$config)$abundance, "double")
  expect_type(LFQData$new(x$data, x$config, setup = TRUE)$data_long()$abundance, "double")
  set_up <- setup_analysis(x$data, x$config)
  set_up$abundance <- as.integer(set_up$abundance)
  expect_type(LFQData$new(set_up, x$config)$data_long()$abundance, "double")
})

test_that("factor and hierarchy columns of class Date or AsIs are used by their labels", {
  withr::local_options(prolfqua.verbose = FALSE)
  x <- types_input(factor_column = "day")
  days <- c(Ctrl = "2024-01-01", A = "2024-02-01", B = "2024-03-01")
  x$data$day <- as.Date(days[x$data$group])
  x$data$proteinID <- I(x$data$proteinID)
  lfq <- LFQData$new(x$data, x$config, setup = TRUE)
  expect_setequal(unique(lfq$data_long()$day_), unname(days))
  expect_type(lfq$data_long()$protein_Id, "character")
  expect_setequal(lfq$get_Stats()$stats()$day_, c(unname(days), "All"))
})

test_that("list columns in configured roles are refused; other list columns are dropped", {
  withr::local_options(prolfqua.verbose = FALSE)
  x <- types_input()
  for (column in c("group", "proteinID", "sample")) {
    data <- tibble::as_tibble(x$data)
    data[[column]] <- as.list(data[[column]])
    expect_error(
      LFQData$new(data, x$config, setup = TRUE),
      paste0("list columns `", column, "`"),
      class = "prolfqua_error_bad_argument"
    )
  }
  data <- tibble::as_tibble(x$data)
  data$notes <- as.list(seq_len(nrow(data)))
  expect_false("notes" %in% names(LFQData$new(data, x$config, setup = TRUE)$data_long()))
})

test_that("data.frame, tibble and data.table input give the same results", {
  skip_if_not_installed("data.table")
  withr::local_options(prolfqua.verbose = FALSE)
  x <- types_input()
  inputs <- list(
    data.frame = x$data,
    tibble = tibble::as_tibble(x$data),
    data.table = data.table::as.data.table(x$data)
  )
  contrasts <- c(AvsCtrl = "group_A - group_Ctrl")
  results <- lapply(inputs, function(data) {
    lfq <- LFQData$new(data, x$config, setup = TRUE)
    expect_s3_class(lfq$data_long(), "tbl_df")
    build_contrast_analysis(lfq, "~ group_", contrasts, method = "lm")$get_contrasts()
  })
  expect_equal(results$tibble, results$data.frame)
  expect_equal(results$data.table, results$data.frame)
  # without setup, the stored table is a tibble too
  set_up <- setup_analysis(x$data, x$config)
  expect_s3_class(LFQData$new(data.table::as.data.table(set_up), x$config)$data_long(), "tbl_df")
})
