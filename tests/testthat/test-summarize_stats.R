# .summarize_stats() against a per-group dplyr::summarize() reference

reference_summarize_stats <- function(lfqdata, factor_key) {
  resp <- lfqdata$response()
  lfqdata$data_long() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(lfqdata$hierarchy_keys(), factor_key)))) |>
    dplyr::summarize(
      nrReplicates = dplyr::n(),
      nrMeasured = sum(!is.na(.data[[resp]])),
      nrNAs = sum(is.na(.data[[resp]])),
      sd = stats::sd(.data[[resp]], na.rm = TRUE),
      var = stats::var(.data[[resp]], na.rm = TRUE),
      meanAbundance = mean(.data[[resp]], na.rm = TRUE),
      medianAbundance = stats::median(.data[[resp]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(factor_key), as.character))
}

stat_columns <- c("nrReplicates", "nrMeasured", "nrNAs", "sd", "var", "meanAbundance", "medianAbundance")

test_that("summarize_stats equals per-group summaries, including empty and single-value groups", {
  sim <- prolfqua::sim_lfq_data_peptide_config(n_proteins = 20, weight_missing = 0.5, seed = 7)
  # one protein without any value, one with a single value per group
  ids <- unique(sim$data$protein_Id)
  sim$data$abundance[sim$data$protein_Id == ids[1]] <- NA
  single <- which(sim$data$protein_Id == ids[2])
  sim$data$abundance[single[duplicated(sim$data$group_[single])]] <- NA
  lfq <- suppressMessages(prolfqua::LFQData$new(sim$data, sim$config))

  for (factor_key in list(lfq$relevant_factor_keys(), NULL)) {
    got <- .summarize_stats(lfq, factor_key = factor_key)
    expected <- reference_summarize_stats(lfq, factor_key)
    keys <- setdiff(colnames(expected), stat_columns)
    expect_identical(as.data.frame(got[keys]), as.data.frame(expected[keys]))
    expect_equal(as.data.frame(got[stat_columns]), as.data.frame(expected[stat_columns]), tolerance = 1e-12)
  }
  empty <- .summarize_stats(lfq)
  empty <- empty[empty$protein_Id == ids[1], ]
  expect_true(all(is.nan(empty$meanAbundance)))
  expect_true(all(is.na(empty$medianAbundance) & is.na(empty$sd)))
  expect_type(empty$nrMeasured, "integer")
})
