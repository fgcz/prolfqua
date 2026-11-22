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
  withr::local_seed(42)
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
  withr::local_seed(42)
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

test_that("compute_contrast matches summary.lm and the hand-computed Wald test next to a non-estimable row", {
  withr::local_seed(3)
  df <- data.frame(
    y = c(rnorm(4, mean = 10), rnorm(4, mean = 12)),
    group = factor(rep(c("A", "B"), each = 4), levels = c("A", "B", "C"))
  )
  m <- lm(y ~ group, data = df)
  # lm drops the empty level C, so groupC is a linfct column the model does not estimate
  linfct <- rbind(B = c(0, 1, 0), C = c(0, 0, 1), mix = c(1, 0.5, 0))
  colnames(linfct) <- c("(Intercept)", "groupB", "groupC")
  res <- compute_contrast(m, linfct, confint = 0.9)

  ref <- summary(m)$coefficients["groupB", ]
  expect_equal(res$estimate[1], unname(ref["Estimate"]))
  expect_equal(res$std.error[1], unname(ref["Std. Error"]))
  expect_equal(res$p.value[1], unname(ref["Pr(>|t|)"]))
  expect_true(all(is.na(res[2, c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")])))

  L <- linfct["mix", 1:2, drop = FALSE]
  est <- unname(drop(L %*% coefficients(m)[1:2]))
  se <- sqrt(drop(L %*% vcov(m)[1:2, 1:2] %*% t(L)))
  q <- stats::qt(0.95, df.residual(m))
  expect_equal(res$estimate[3], est)
  expect_equal(res$std.error[3], se)
  expect_equal(c(res$conf.low[3], res$conf.high[3]), est + c(-1, 1) * q * se)
  expect_equal(res$df, rep(df.residual(m), 3))
  expect_equal(res$sigma, rep(sigma(m), 3))
})

test_that("linfct_matrix_contrasts names an unnamed contrast contrast_<i>", {
  linfct <- linfct_from_model(sim_make_model_lm("parallel3"), as_list = FALSE)
  res <- expect_no_warning(linfct_matrix_contrasts(linfct, c("group_A - group_Ctrl", AB = "group_A - group_B")))
  expect_identical(rownames(res), c("contrast_1", "AB"))
  expect_equal(res["contrast_1", ], linfct["group_A", ] - linfct["group_Ctrl", ])
})

test_that("an infinite variance of an unused coefficient does not reach the contrast (rfit on small groups)", {
  linfct <- rbind(uses_b = c(0, 1), uses_a = c(1, 0))
  V <- matrix(c(Inf, 0, 0, 2), 2)
  expect_equal(prolfqua:::.contrast_variance(linfct, V), c(2, Inf))
})
