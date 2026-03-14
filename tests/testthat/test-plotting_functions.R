test_that("plot_hierarchies_boxplot", {
  istar <- sim_lfq_data_peptide_config()
  analysis <- istar$data
  config <- istar$config
  config$hierarchyDepth
  config$hierarchy_keys_depth()
  #'
  xnested <- analysis |>
    dplyr::group_by_at(config$hierarchy_keys_depth()) |>
    tidyr::nest()
  #'
  p <- plot_hierarchies_boxplot(
    xnested$data[[3]],
    xnested$protein_Id[[3]],
    config,
    facet_grid_on = tail(config$hierarchy_keys(), 1)
  )
  expect_true("ggplot" %in% class(p))

  p <- plot_hierarchies_boxplot(xnested$data[[3]], xnested$protein_Id[[3]], config)
  expect_true("ggplot" %in% class(p))
  p <- plot_hierarchies_boxplot(
    xnested$data[[3]],
    xnested$protein_Id[[3]],
    config,
    beeswarm = FALSE
  )
  expect_true("ggplot" %in% class(p))
})
