# Per-protein QC figures: plots that do not carry the whole data set, writers that render one plot at a time

log2_peptide_lfqdata <- function(n_proteins) {
  sim <- prolfqua::sim_lfq_data_peptide_config(n_proteins = n_proteins, seed = 5)
  lfq <- suppressMessages(prolfqua::LFQData$new(sim$data, sim$config))
  suppressMessages(lfq$get_Transformer()$log2()$lfq)
}

test_that("one aggregation plot and one boxplot do not grow with the number of proteins", {
  plot_sizes <- function(n_proteins) {
    lfq <- log2_peptide_lfqdata(n_proteins)
    agg <- estimate_intensity(lfq, method = "medpolish")
    lfq_agg <- suppressMessages(prolfqua::LFQData$new(agg$data, agg$config))
    c(
      estimate = length(serialize(plot_estimate(lfq, lfq_agg)$plots[[1]], NULL)),
      boxplot = length(serialize(lfq$get_Plotter()$boxplots()$boxplot[[1]], NULL))
    )
  }
  small <- plot_sizes(4)
  large <- plot_sizes(60)
  # before the fix one aggregation plot of 40 proteins serialised to 268 MB
  expect_true(all(large < 1.5 * small))
})

test_that("write_plots renders one plot per protein and passes show.legend", {
  lfq <- log2_peptide_lfqdata(8)
  agg <- suppressMessages(prolfqua::AggregateMedpolish$new(lfq, "protein"))
  suppressMessages(agg$aggregate())
  rendered <- list()
  local_mocked_bindings(.render_plot_to_device = function(fig) rendered[[length(rendered) + 1]] <<- fig)
  path <- agg$write_plots(withr::local_tempdir(), show.legend = TRUE)
  expect_true(file.exists(path))
  expect_length(rendered, length(unique(agg$lfq_agg$data_long()$protein_Id)))
  feature_layers <- rendered[[1]]$layers[1:2]
  expect_true(all(vapply(feature_layers, function(layer) isTRUE(layer$show.legend), logical(1))))
})

test_that("write_boxplots renders one plot per protein", {
  lfq <- log2_peptide_lfqdata(8)
  rendered <- 0L
  local_mocked_bindings(.render_plot_to_device = function(fig) rendered <<- rendered + 1L)
  path <- suppressMessages(lfq$get_Plotter()$write_boxplots(withr::local_tempdir()))
  expect_true(file.exists(path))
  expect_identical(rendered, length(unique(lfq$data_long()$protein_Id)))
})
