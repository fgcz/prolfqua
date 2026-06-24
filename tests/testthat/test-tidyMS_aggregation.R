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

# Regression tests for .rlm_estimate (TODO_prolfqua_prolfquapp_code_review.md).
# The review flagged "three different column schemas". The column SET is in fact identical in all
# three branches and the per-sample key is never dropped (consumers bind by name), so there were no
# wrong results -- but column ORDER diverged and only the main branch re-established all samples.
# These pin one consistent schema + full sample coverage across every branch.

test_that(".rlm_estimate returns one consistent column schema across all three branches", {
  expected <- c("samples", "mean.response", "lmrob", "weights")

  # branch 1: single feature, multiple samples
  df_feat1 <- data.frame(samples = letters[1:4], feature = "F1", response = c(10, 11, 12, 13))
  res_feat1 <- prolfqua:::.rlm_estimate(df_feat1, "response", "feature", "samples")

  # branch 2: multiple features, single sample
  df_samp1 <- data.frame(samples = "a", feature = c("F1", "F2", "F3"), response = c(10, 11, 12))
  res_samp1 <- prolfqua:::.rlm_estimate(df_samp1, "response", "feature", "samples")

  # main path: multiple features and samples (mirrors the rlm_estimate roxygen example)
  set.seed(42)
  xx <- data.frame(response = rnorm(20, 0, 10), feature = rep(LETTERS[1:5], 4), samples = rep(letters[1:4], 5))
  res_main <- prolfqua:::.rlm_estimate(xx, "response", "feature", "samples")

  # same column SET in every branch (always true -- documents the no-wrong-results finding)
  expect_setequal(colnames(res_feat1), expected)
  expect_setequal(colnames(res_samp1), expected)
  expect_setequal(colnames(res_main), expected)

  # identical column ORDER in every branch (the defect: order diverged before the fix)
  expect_identical(colnames(res_feat1), expected)
  expect_identical(colnames(res_samp1), expected)
  expect_identical(colnames(res_main), expected)
})

test_that(".rlm_estimate re-establishes every input sample in all branches", {
  # single-feature unit with explicit NA-response rows -- samples b and d must NOT be dropped
  df_na <- data.frame(samples = letters[1:4], feature = "F1", response = c(10, NA, 12, NA))
  res_na <- prolfqua:::.rlm_estimate(df_na, "response", "feature", "samples")
  expect_setequal(res_na$samples, letters[1:4])
  expect_equal(nrow(res_na), 4L)

  # main path already re-establishes all samples -- pin that it stays true
  set.seed(7)
  xx <- data.frame(response = rnorm(20, 0, 10), feature = rep(LETTERS[1:5], 4), samples = rep(letters[1:4], 5))
  res_main <- prolfqua:::.rlm_estimate(xx, "response", "feature", "samples")
  expect_setequal(res_main$samples, letters[1:4])
})

test_that("rlm and medpolish aggregation produce the same (protein, sample) rows (boundary invariant)", {
  dd <- prolfqua::sim_lfq_data_peptide_config()
  lfq <- prolfqua::LFQData$new(dd$data, dd$config)
  lfq <- lfq$get_Transformer()$log2()$lfq

  bbMed <- suppressMessages(prolfqua::estimate_intensity(lfq, .func = prolfqua::medpolish_estimate_dfconfig))
  bbRob <- suppressMessages(prolfqua::estimate_intensity(lfq, .func = prolfqua::rlm_estimate_dfconfig))

  keycols <- c(bbMed$config$hierarchy_keys(), bbMed$config$sample_name)
  keyM <- dplyr::distinct(bbMed$data, dplyr::across(dplyr::all_of(keycols)))
  keyR <- dplyr::distinct(bbRob$data, dplyr::across(dplyr::all_of(keycols)))
  expect_equal(nrow(dplyr::anti_join(keyM, keyR, by = keycols)), 0L)
  expect_equal(nrow(dplyr::anti_join(keyR, keyM, by = keycols)), 0L)

  # rlm carries its estimate, fitting weights, and the nr_children count
  expect_true(all(c("lmrob", "weights") %in% colnames(bbRob$data)))
  expect_true(bbRob$config$nr_children %in% colnames(bbRob$data))
})
