# wide group-mean table consumed by get_contrast (one row per protein, one column per group)
.contrast_wide_fixture <- function() {
  tibble::tibble(
    protein_Id = c("p1", "p2"),
    group_A = c(10, 20),
    group_B = c(4, 5),
    group_Ctrl = c(1, 2)
  )
}

test_that("get_contrast preserves contrast-side order for a simple contrast", {
  data <- .contrast_wide_fixture()
  res <- suppressMessages(prolfqua::get_contrast(data, "protein_Id", c(AvsB = "group_A - group_B")))

  expect_equal(res$group_1, c(10, 20))
  expect_equal(res$group_2, c(4, 5))
  expect_equal(res$estimate, c(6, 15))
  expect_equal(unique(res$group_1_name), "group_A")
  expect_equal(unique(res$group_2_name), "group_B")
})

test_that("get_contrast errors loudly when a contrast level is absent", {
  data <- .contrast_wide_fixture()
  expect_error(
    suppressMessages(prolfqua::get_contrast(data, "protein_Id", c(bad = "group_A - group_Xyz")))
  )
})

test_that("get_contrast rejects a contrast that is not a difference 'LHS - RHS'", {
  data <- .contrast_wide_fixture()
  expect_error(
    suppressMessages(prolfqua::get_contrast(data, "protein_Id", c(bad = "group_A + group_B"))),
    "LHS - RHS"
  )
})

test_that("get_contrast derives group_1/group_2 from the contrast sides for averaging contrasts", {
  data <- .contrast_wide_fixture()
  res <- suppressMessages(
    prolfqua::get_contrast(data, "protein_Id", c(AvB_vs_Ctrl = "(group_A + group_B)/2 - group_Ctrl"))
  )

  expect_equal(res$group_1, c((10 + 4) / 2, (20 + 5) / 2))
  expect_equal(res$group_2, c(1, 2))
  expect_equal(res$estimate, c((10 + 4) / 2 - 1, (20 + 5) / 2 - 2))
})

test_that("get_contrast resolves nested contrasts that reference earlier contrast names", {
  # Regression: a later contrast may reference the names of earlier contrasts as
  # groups (e.g. interaction-style "T_C_gv_WT - T_C_gv_KO"). Each contrast is
  # materialized under its own name so subsequent contrasts can resolve it.
  data <- .contrast_wide_fixture()
  res <- suppressMessages(prolfqua::get_contrast(
    data,
    "protein_Id",
    c(
      AvsB = "group_A - group_B",
      AvsCtrl = "group_A - group_Ctrl",
      nested = "AvsB - AvsCtrl"
    )
  ))

  nested <- dplyr::filter(res, .data[["contrast"]] == "nested")
  avsb <- c(10 - 4, 20 - 5)
  avsctrl <- c(10 - 1, 20 - 2)
  expect_equal(nested$group_1, avsb)
  expect_equal(nested$group_2, avsctrl)
  expect_equal(nested$estimate, avsb - avsctrl)
  expect_equal(unique(nested$group_1_name), "AvsB")
  expect_equal(unique(nested$group_2_name), "AvsCtrl")
})
