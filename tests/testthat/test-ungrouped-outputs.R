# Exported functions return ungrouped tibbles: a grouped result makes later user mutate()/summarise() calls run
# per group without asking for it.

.ungrouped_lfq <- function() {
  istar <- sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE, seed = 1)
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
  by_group <- summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), factors = "group_")
  expect_false(dplyr::is_grouped_df(by_group))
  counts <- HierarchyCountsSample$new(lfq$data_long(), lfq$get_config())
  expect_false(dplyr::is_grouped_df(counts$wide()))
  expect_false(dplyr::is_grouped_df(counts$long()))
})

test_that("missingness_per_condition() and _cumsum() return ungrouped data", {
  lfq <- .ungrouped_lfq()
  expect_false(dplyr::is_grouped_df(missingness_per_condition(lfq)$data))
  expect_false(dplyr::is_grouped_df(missingness_per_condition_cumsum(lfq)$data))
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
