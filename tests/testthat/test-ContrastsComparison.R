# ContrastsComparison: aligning several methods on their subject and contrast pairs

comparison_tables <- function() {
  # P1 up in a, down in b; P2 up in a, not significant in b; P3 not significant in a, up in b; P4 without effect in a,
  # absent from b; P5 up in a, absent from b
  a <- data.frame(
    protein_Id = c("P1", "P2", "P3", "P4", "P5"),
    contrast = "c1",
    diff = c(1, 1, 0.2, NA, 2),
    statistic = c(5, 5, 0.5, NA, 6),
    p.value = c(0.001, 0.001, 0.6, NA, 0.001),
    FDR = c(0.01, 0.01, 0.7, NA, 0.01)
  )
  b <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    contrast = "c1",
    diff = c(-1, 0.5, 1),
    statistic = c(-5, 1, 5),
    p.value = c(0.001, 0.4, 0.001),
    FDR = c(0.01, 0.5, 0.01)
  )
  list(a = ContrastsTable$new(a, subject_id = "protein_Id"), b = ContrastsTable$new(b, subject_id = "protein_Id"))
}

test_that("ContrastsComparison aligns the methods on every reported pair with a status", {
  comparison <- ContrastsComparison$new(comparison_tables())
  aligned <- comparison$aligned()
  expect_equal(nrow(aligned), 10)
  status <- stats::setNames(aligned$status, paste(aligned$protein_Id, aligned$method))
  expect_equal(
    unname(status[c("P1 a", "P1 b", "P2 b", "P3 a", "P4 a", "P5 b")]),
    c("up", "down", "not significant", "not significant", "not estimated", "not estimated")
  )
})

test_that("ContrastsComparison flags sign, significance and coverage conflicts", {
  conflicts <- ContrastsComparison$new(comparison_tables())$conflicts()
  flags <- conflicts[order(conflicts$protein_Id), c("sign_conflict", "significance_conflict", "coverage_conflict")]
  expect_equal(flags$sign_conflict, c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(flags$significance_conflict, c(FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(flags$coverage_conflict, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(conflicts$not_estimated[conflicts$protein_Id == "P4"], 2)
})

test_that("ContrastsComparison counts the pairs per set of estimating methods", {
  coverage <- ContrastsComparison$new(comparison_tables())$coverage()
  expect_equal(stats::setNames(coverage$pairs, coverage$estimated_by), c("a & b" = 3, a = 1, none = 1))
})

test_that("ContrastsComparison reads other column names through the ContrastConfiguration", {
  config <- ContrastConfiguration$new(
    subject_id = "protein_id",
    effect_col = "log_fc",
    score_col = "t_statistic",
    pvalue_col = "p_value",
    fdr_col = "p_value_adjusted"
  )
  tool <- function(log_fc, p_value) {
    table <- data.frame(
      protein_id = c("P1", "P2", "P3"),
      contrast = "c1",
      log_fc = log_fc,
      t_statistic = log_fc * 3,
      p_value = p_value,
      p_value_adjusted = p_value * 2
    )
    ContrastsTable$new(table, subject_id = "protein_id", config = config)
  }
  comparison <- ContrastsComparison$new(list(
    x = tool(c(1, -1, 0.5), c(0.001, NA, 0.2)),
    y = tool(c(1, -1, NA), c(0.001, 0.001, NA))
  ))
  aligned <- comparison$aligned()
  expect_equal(aligned$effect[aligned$method == "x"], c(1, -1, 0.5))
  expect_equal(aligned$FDR[aligned$method == "y"], c(0.002, 0.002, NA))
  expect_true(all(is.na(aligned$estimate_type)))

  # P2 has a fold change but no p-value in x
  expect_equal(comparison$pairs("x", columns = "effect")$protein_id, c("P1", "P2", "P3"))
  expect_equal(comparison$pairs("x", columns = c("effect", "p.value"))$protein_id, c("P1", "P3"))
  expect_equal(comparison$pairs("y", columns = "statistic")$protein_id, c("P1", "P2"))
})

test_that("ContrastsComparison finds the pairs an imputing method adds to those of lm", {
  sim <- sim_lfq_data_protein_config(n_proteins = 30, weight_missing = 0, seed = 1)
  protein <- sort(unique(sim$data$protein_Id))[1]
  data <- sim$data[!(sim$data$protein_Id == protein & sim$data$group_ == "A"), ]
  lfqdata <- LFQData$new(data, sim$config)
  contrasts <- c(A_vs_Ctrl = "group_A - group_Ctrl", B_vs_Ctrl = "group_B - group_Ctrl")
  fits <- lapply(c(lm = "lm", limma_impute = "limma_impute"), function(method) {
    suppressMessages(build_contrast_analysis(lfqdata, "~ group_", contrasts, method = method))
  })
  comparison <- suppressWarnings(ContrastsComparison$new(fits))

  coverage <- comparison$coverage()
  expect_equal(coverage$pairs[coverage$estimated_by == "limma_impute"], 1)
  lm_pairs <- comparison$pairs("lm")
  expect_false(any(lm_pairs$protein_Id == protein & lm_pairs$contrast == "A_vs_Ctrl"))
  added <- comparison$aligned()
  added <- added[added$method == "limma_impute" & added$protein_Id == protein & added$contrast == "A_vs_Ctrl", ]
  expect_equal(added$estimate_type, "lod_imputed")
})

test_that("ContrastsComparison refuses methods that do not share subjects and contrasts", {
  tables <- comparison_tables()
  expect_error(ContrastsComparison$new(unname(tables)), "unique names")
  expect_error(ContrastsComparison$new(tables["a"]), "at least two")
  expect_error(ContrastsComparison$new(tables, fdr_threshold = 5), "between 0 and 1")

  other_subject <- tables$b$get_contrasts()
  names(other_subject)[1] <- "protein_id"
  tables$b <- ContrastsTable$new(other_subject, subject_id = "protein_id")
  expect_error(ContrastsComparison$new(tables), "share their subject columns")

  tables <- comparison_tables()
  other_contrast <- tables$b$get_contrasts()
  other_contrast$contrast <- "c2"
  tables$b <- ContrastsTable$new(other_contrast, subject_id = "protein_Id")
  expect_error(ContrastsComparison$new(tables), "same contrast names")
})

test_that("ContrastsComparison draws an UpSet plot of the estimated pairs", {
  skip_if_not_installed("UpSetR")
  # UpSetR still calls ggplot2 functions deprecated since ggplot2 3.0 and 3.4
  withr::local_options(lifecycle_verbosity = "quiet")
  expect_s3_class(ContrastsComparison$new(comparison_tables())$upset(), "upset")
})
