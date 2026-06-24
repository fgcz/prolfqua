# Regression tests for nr_obs_sample / nr_children_experiment.
#
# These pin the contract that a code review (TODO_prolfqua_prolfquapp_code_review.md) misread as a bug:
# nr_children_experiment maxes the per-SAMPLE SUMMED child count, not the un-aggregated source column.
# The misread came from nr_obs_sample's default `new_child = nr_children_col` reusing the source name.

.nr_child_fixture <- function() {
  dd <- prolfqua::sim_lfq_data_peptide_config()
  prolfqua::LFQData$new(dd$data, dd$config)
}

test_that("nr_obs_sample sums nr_children per (hierarchy, sample) into the new_child column", {
  lfq <- .nr_child_fixture()
  d <- lfq$data_long()
  nr_col <- lfq$nr_children_col()

  out <- prolfqua::nr_obs_sample(
    d,
    lfq$response(),
    lfq$relevant_hierarchy_keys(),
    lfq$file_name(),
    nr_children_col = nr_col,
    new_child = "summed_children"
  )

  # output carries the requested column name (not the source name)
  expect_true("summed_children" %in% colnames(out))
  expect_false(nr_col %in% setdiff(colnames(out), c(lfq$relevant_hierarchy_keys(), lfq$file_name())))

  # values equal an independent per-(protein, sample) sum of nr_children over observed rows
  ref <- d[!is.na(d[[lfq$response()]]), ] |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(lfq$relevant_hierarchy_keys(), lfq$file_name())))) |>
    dplyr::summarize(ref = sum(.data[[nr_col]], na.rm = TRUE), .groups = "drop")
  chk <- dplyr::inner_join(out, ref, by = c(lfq$relevant_hierarchy_keys(), lfq$file_name()))
  expect_equal(chk$summed_children, chk$ref)
})

test_that("nr_children_experiment maxes the per-sample SUMMED count, not the un-aggregated column", {
  lfq <- .nr_child_fixture()
  d <- lfq$data_long()
  nr_col <- lfq$nr_children_col()

  # In peptide-level sim data every observed per-row nr_children == 1.
  observed <- d[!is.na(d[[lfq$response()]]), ]
  expect_equal(unique(observed[[nr_col]]), 1)

  xd <- prolfqua::nr_children_experiment(
    d,
    lfq$response(),
    lfq$relevant_hierarchy_keys(),
    lfq$file_name(),
    nr_children_col = nr_col
  )

  # Killer assertion: if the function maxed the un-aggregated all-1 column, every value would be 1.
  # It is > 1 only because it maxes the per-sample SUM (peptides observed per protein in a sample).
  expect_gt(max(xd$nr_child_exp), 1)
  expect_equal(min(xd$nr_child_exp), 1) # roxygen @examples invariant

  # Cross-check against an independent max_over_samples(sum_per_sample(nr_children)).
  ref <- d[!is.na(d[[lfq$response()]]), ] |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(lfq$relevant_hierarchy_keys(), lfq$file_name())))) |>
    dplyr::summarize(s = sum(.data[[nr_col]], na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(lfq$relevant_hierarchy_keys()))) |>
    dplyr::summarize(ref = max(s), .groups = "drop")
  chk <- dplyr::inner_join(xd, ref, by = lfq$relevant_hierarchy_keys())
  expect_equal(chk$nr_child_exp, chk$ref)
})
