test_that("plot_hierarchies_boxplot", {
  istar <- sim_lfq_data_peptide_config()
  lfq <- LFQData$new(istar$data, istar$config)

  xnested <- lfq$data_long() |>
    dplyr::group_by_at(lfq$relevant_hierarchy_keys()) |>
    tidyr::nest()

  p <- plot_hierarchies_boxplot(
    xnested$data[[3]],
    xnested$protein_Id[[3]],
    lfq,
    facet_grid_on = tail(lfq$hierarchy_keys(), 1)
  )
  expect_true("ggplot" %in% class(p))

  p <- plot_hierarchies_boxplot(xnested$data[[3]], xnested$protein_Id[[3]], lfq)
  expect_true("ggplot" %in% class(p))
  p <- plot_hierarchies_boxplot(
    xnested$data[[3]],
    xnested$protein_Id[[3]],
    lfq,
    beeswarm = FALSE
  )
  expect_true("ggplot" %in% class(p))
})
