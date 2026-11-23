# Exported functions return ungrouped tibbles: a grouped result makes later user mutate()/summarise() calls run
# per group without asking for it.

.ungrouped_lfq <- function() {
  istar <- sim_lfq_data_peptide_config(n_proteins = 10, with_missing = TRUE, seed = 1)
  LFQData$new(istar$data, istar$config)
}

test_that("adjust_p_values() returns ungrouped data with per-contrast BH adjustment", {
  df <- data.frame(contrast = rep(c("a", "b"), each = 5), p.value = c(1:5, 6:10) / 20)
  res <- adjust_p_values(df)
  expect_false(dplyr::is_grouped_df(res))
  expect_equal(res$FDR, stats::ave(df$p.value, df$contrast, FUN = function(p) stats::p.adjust(p, "BH")))
})

test_that("summarize_hierarchy() and HierarchyCountsSample return ungrouped data", {
  lfq <- .ungrouped_lfq()
  expect_false(dplyr::is_grouped_df(lfq$summarize_hierarchy()))
  by_group <- .summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), factors = "group_")
  expect_false(dplyr::is_grouped_df(by_group))
  counts <- HierarchyCountsSample$new(lfq$data_long(), lfq$get_config())
  expect_false(dplyr::is_grouped_df(counts$wide()))
  expect_false(dplyr::is_grouped_df(counts$long()))
})

test_that("missingness_per_condition() and _cumsum() return ungrouped data", {
  lfq <- .ungrouped_lfq()
  expect_false(dplyr::is_grouped_df(.missingness_per_condition(lfq)$data))
  expect_false(dplyr::is_grouped_df(.missingness_per_condition_cumsum(lfq)$data))
})

test_that("barplot_threshold() summaries and duplicate counts are ungrouped", {
  lfq <- .ungrouped_lfq()
  lfq <- lfq$get_Transformer()$log2()$lfq
  agg <- lfq$get_Aggregator()
  agg$aggregate()
  facade <- build_contrast_analysis(
    agg$lfq_agg,
    "~ group_",
    c(A_vs_Ctrl = "group_A - group_Ctrl", B_vs_Ctrl = "group_B - group_Ctrl"),
    method = "lm"
  )
  cp <- ContrastsPlotter$new(facade$get_contrasts(), "protein_Id", volcano = list(list(score = "FDR", thresh = 1)))
  expect_false(dplyr::is_grouped_df(cp$barplot_threshold()$FDR$summary))
  expect_false(dplyr::is_grouped_df(.setup_duplicate_counts(lfq$data_long(), lfq$get_config())))
})

test_that("table_factors_size returns an ungrouped tibble for two factors", {
  sim <- sim_lfq_data_2factor_config(n_proteins = 5)
  lfq <- LFQData$new(sim$data, sim$config)
  res <- .table_factors_size(
    lfq$data_long(),
    lfq$file_name(),
    lfq$sample_name(),
    lfq$factor_keys(),
    lfq$relevant_factor_keys()
  )
  expect_false(dplyr::is_grouped_df(res))
  expect_identical(sum(res$n), dplyr::n_distinct(lfq$data_long()[[lfq$sample_name()]]))
})

test_that("summarize_stats_quantiles returns ungrouped quantiles per factor level (review D4)", {
  sim <- sim_lfq_data_2factor_config(n_proteins = 20)
  lfq <- LFQData$new(sim$data, sim$config)
  stats <- .summarize_stats(lfq)
  probs <- c(0.25, 0.5, 0.9)
  res <- .summarize_stats_quantiles(stats, lfq$relevant_factor_keys(), probs = probs)
  expect_false(dplyr::is_grouped_df(res$long))
  expect_identical(colnames(res$long), c(lfq$relevant_factor_keys(), "probs", "quantiles"))
  levels_n <- nrow(dplyr::distinct(stats, dplyr::across(dplyr::all_of(lfq$relevant_factor_keys()))))
  expect_identical(nrow(res$long), levels_n * length(probs))
  one <- dplyr::filter(stats, .data$Treatment == "A", .data$Background == "X", !is.na(.data$sd))
  expected <- stats::quantile(one$sd, probs, na.rm = TRUE)
  got <- dplyr::filter(res$long, .data$Treatment == "A", .data$Background == "X")
  expect_equal(unname(got$quantiles), unname(expected))
  expect_identical(colnames(res$wide), c("probs", "A_X", "A_Z", "B_X", "B_Z"))
})
