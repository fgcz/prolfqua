# null_split_check()

.ns_lfq <- function(Nprot = 30, with_missing = FALSE) {
  bb <- sim_lfq_data_protein_config(Nprot = Nprot, with_missing = with_missing)
  lfq <- LFQData$new(bb$data, bb$config)
  lfq$is_transformed(TRUE)
  lfq
}

test_that("one row per split with summary and list-columns", {
  res <- null_split_check(.ns_lfq(), method = "limma", n_splits = 4, seed = 1)
  expect_equal(nrow(res), 4)
  expect_named(
    res,
    c("split", "n", "p_below_0.05", "p_below_0.01", "n_FDR_below_0.05", "ks_D", "assignment", "contrasts")
  )
  expect_true(all(res$p_below_0.05 >= 0 & res$p_below_0.05 <= 1))
  expect_true(all(res$ks_D >= 0 & res$ks_D <= 1))
  ct <- res$contrasts[[1]]
  expect_true(all(c("protein_Id", "diff", "p.value", "FDR") %in% colnames(ct)))
  expect_equal(res$n[1], sum(!is.na(ct$p.value)))
  expect_equal(res$p_below_0.05[1], mean(ct$p.value < 0.05, na.rm = TRUE))
})

test_that("splits are balanced within strata and distinct up to swapping A and B", {
  lfq <- .ns_lfq()
  res <- null_split_check(lfq, method = "limma", n_splits = 6, seed = 2)
  groups <- dplyr::distinct(lfq$data_long(), sampleName, group_)
  for (a in res$assignment) {
    counts <- table(dplyr::inner_join(a, groups, by = "sampleName")[, c("group_", "pseudo_group_")])
    expect_true(all(counts == 2))
  }
  canonical <- vapply(
    res$assignment,
    function(a) {
      g <- a$pseudo_group_
      if (g[1] == "B") {
        g <- ifelse(g == "A", "B", "A")
      }
      paste(g, collapse = "")
    },
    character(1)
  )
  expect_equal(anyDuplicated(canonical), 0L)
})

test_that("one group of four samples has three distinct splits", {
  lfq <- .ns_lfq()
  lfq$set_data(lfq$data_long() |> dplyr::filter(group_ == "Ctrl"))
  expect_message(
    res <- null_split_check(lfq, method = "limma", n_splits = 10, strata = character(0), seed = 3),
    "only 3 distinct splits"
  )
  expect_equal(nrow(res), 3)
})

test_that("same seed gives the same result", {
  lfq <- .ns_lfq()
  r1 <- null_split_check(lfq, method = "limma", n_splits = 3, n_features = 10, seed = 5)
  r2 <- null_split_check(lfq, method = "limma", n_splits = 3, n_features = 10, seed = 5)
  expect_identical(r1, r2)
})

test_that("n_features restricts the tested proteins", {
  res <- null_split_check(.ns_lfq(), method = "limma", n_splits = 2, n_features = 10, seed = 6)
  for (ct in res$contrasts) {
    expect_lte(length(unique(ct$protein_Id)), 10)
  }
  ids <- lapply(res$contrasts, function(ct) sort(unique(ct$protein_Id)))
  expect_identical(ids[[1]], ids[[2]])
})

test_that("invalid strata and splits raise", {
  lfq <- .ns_lfq()
  expect_error(null_split_check(lfq, strata = "batch"), class = "prolfqua_error_bad_argument")
  expect_error(null_split_check(lfq, n_splits = 0), class = "prolfqua_error_bad_argument")
  one <- lfq$get_copy()
  d <- one$data_long()
  one$set_data(d[!(d$group_ == "A" & d$sampleName != "A_V1"), ])
  expect_error(null_split_check(one, method = "limma"), "at least two samples")
})

test_that("p-values of a no-effect simulation are roughly calibrated", {
  lfq <- .ns_lfq(Nprot = 200)
  res <- null_split_check(lfq, method = "limma", n_splits = 10, seed = 7)
  expect_lt(stats::median(res$p_below_0.05), 0.15)
  expect_lt(stats::median(res$n_FDR_below_0.05), 5)
})
