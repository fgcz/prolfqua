test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("plot_hierarchies_boxplot",{

  istar <- sim_lfq_data_peptide_config()
  analysis <- istar$data
  config <- istar$config
  config$table$hierarchyDepth
  config$table$hierarchy_keys_depth()
  #'
  xnested <- analysis |>
    dplyr::group_by_at(config$table$hierarchy_keys_depth()) |> tidyr::nest()
  #'
  p <- plot_hierarchies_boxplot(xnested$data[[3]],
                                xnested$protein_Id[[3]],
                                config,
                                facet_grid_on =  tail(config$table$hierarchy_keys(),1))
  stopifnot("ggplot" %in% class(p))

  p <- plot_hierarchies_boxplot(xnested$data[[3]],
                                xnested$protein_Id[[3]],
                                config )
  stopifnot("ggplot" %in% class(p))
  p <- plot_hierarchies_boxplot(
    xnested$data[[3]],
    xnested$protein_Id[[3]],
    config,
    beeswarm = FALSE )
  stopifnot("ggplot" %in% class(p))

}

)
