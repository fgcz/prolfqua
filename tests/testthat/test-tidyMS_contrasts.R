# linfct_matrix_contrasts and compute_contrast (per-protein Wald contrast helpers)

test_that("linfct_matrix_contrasts resolves self-referencing contrasts", {
  Contr <- c(
    "IntoflintoA" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
    "IntoflintoB" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`",
    "IntoflintoX" = "`TreatmentA:BackgroundX` - `TreatmentB:BackgroundX`",
    "IntoflintoZ" = "`TreatmentA:BackgroundZ` - `TreatmentB:BackgroundZ`",
    "interactXZ" = "IntoflintoX - IntoflintoZ",
    "interactAB" = "IntoflintoA - IntoflintoB"
  )
  set.seed(42)
  additive <- linfct_matrix_contrasts(linfct_from_model(sim_make_model_lm("factors"), as_list = FALSE), Contr)
  expect_equal(sum(additive["interactXZ", ]), 0)
  expect_equal(sum(additive["interactAB", ]), 0)
  interaction <- linfct_matrix_contrasts(linfct_from_model(sim_make_model_lm("interaction"), as_list = FALSE), Contr)
  expect_equal(sum(interaction["interactXZ", ]), 1)
  expect_equal(sum(interaction["interactAB", ]), 1)
})

test_that("linfct_matrix_contrasts reports failed contrasts and keeps the valid ones", {
  linfct <- linfct_from_model(sim_make_model_lm("parallel3"), as_list = FALSE)
  contrasts <- c(valid = "group_A - group_Ctrl", invalid = "nonexistent_group - group_Ctrl")
  expect_warning(
    expect_message(
      res <- linfct_matrix_contrasts(linfct, contrasts, p.message = TRUE),
      "valid=group_A - group_Ctrl"
    ),
    "computed 1/2 contrasts; failed 1: invalid"
  )
  expect_identical(rownames(res), "valid")
})

test_that("compute_contrast returns NA for contrasts on non-estimable coefficients", {
  set.seed(42)
  df <- data.frame(
    y = c(rnorm(3, mean = 10), rnorm(3, mean = 12)),
    group = factor(c("A", "A", "A", "B", "B", "B"), levels = c("A", "B", "C"))
  )
  m <- lm(y ~ group, data = df)
  linfct <- matrix(0, nrow = 2, ncol = 3, dimnames = list(c("AvsB", "AvsC"), c("(Intercept)", "groupB", "groupC")))
  linfct["AvsB", ] <- c(0, -1, 0)
  linfct["AvsC", ] <- c(0, 0, -1)
  res <- compute_contrast(m, linfct, confint = 0.95)
  expect_identical(res$lhs, c("AvsB", "AvsC"))
  expect_false(anyNA(res[1, c("estimate", "std.error", "p.value")]))
  expect_true(all(is.na(res[2, c("estimate", "std.error", "p.value")])))

  # B is aliased with A, so Bb2/Bb3 are NA; +1/-1 weights on them cancel but stay non-estimable.
  d <- data.frame(
    y = c(1.2, -0.3, 0.8, 0.1, 2.1, -1.0, 0.4, 0.9, -0.6, 1.5, 0.2, -0.4),
    A = factor(rep(c("a1", "a2", "a3"), each = 4)),
    B = factor(rep(c("b1", "b2", "b3"), each = 4))
  )
  m <- lm(y ~ A + B, data = d)
  na_names <- names(which(is.na(coefficients(m))))
  expect_gte(length(na_names), 2)
  linfct <- matrix(0, nrow = 1, ncol = length(coefficients(m)), dimnames = list("cancel", names(coefficients(m))))
  linfct[1, na_names[1:2]] <- c(1, -1)
  res <- compute_contrast(m, linfct, confint = 0.95)
  expect_true(is.na(res$estimate))
  expect_true(is.na(res$p.value))
})
