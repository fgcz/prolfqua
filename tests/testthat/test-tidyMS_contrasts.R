# .linfct_matrix_contrasts and .compute_estimable_contrast (per-protein Wald contrast helpers)

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
  additive <- .linfct_matrix_contrasts(.linfct_from_model(.sim_make_model_lm("factors"), as_list = FALSE), Contr)
  expect_equal(sum(additive["interactXZ", ]), 0)
  expect_equal(sum(additive["interactAB", ]), 0)
  interaction <- .linfct_matrix_contrasts(.linfct_from_model(.sim_make_model_lm("interaction"), as_list = FALSE), Contr)
  expect_equal(sum(interaction["interactXZ", ]), 1)
  expect_equal(sum(interaction["interactAB", ]), 1)
})

test_that("linfct_matrix_contrasts reports failed contrasts and keeps the valid ones", {
  linfct <- .linfct_from_model(.sim_make_model_lm("parallel3"), as_list = FALSE)
  contrasts <- c(valid = "group_A - group_Ctrl", invalid = "nonexistent_group - group_Ctrl")
  expect_warning(
    expect_message(
      res <- .linfct_matrix_contrasts(linfct, contrasts, p.message = TRUE),
      "valid=group_A - group_Ctrl"
    ),
    "computed 1/2 contrasts; failed 1: invalid"
  )
  expect_identical(rownames(res), "valid")
})

test_that("linfct_matrix_contrasts refuses calls and scalars and drops unparsable or non-numeric contrasts", {
  linfct <- .linfct_from_model(.sim_make_model_lm("parallel3"), as_list = FALSE)
  # a constant would be recycled over the coefficients and shift the rows after it
  expect_error(
    .linfct_matrix_contrasts(linfct, c(one = "1")),
    "at least one group",
    class = "prolfqua_error_bad_argument"
  )
  expect_error(
    .linfct_matrix_contrasts(linfct, c(call = "group_A - mean(group_Ctrl)")),
    "`mean`",
    class = "prolfqua_error_bad_argument"
  )
  expect_warning(
    res <- .linfct_matrix_contrasts(linfct, c(valid = "group_A - group_Ctrl", broken = "group_A -")),
    "failed 1: broken"
  )
  expect_identical(rownames(res), "valid")
  expect_warning(
    res <- .linfct_matrix_contrasts(linfct, c(text = "group_A - 'Ctrl'", valid = "group_A - group_Ctrl")),
    "non-numeric argument"
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
  res <- .compute_estimable_contrast(m, linfct, confint = 0.95)
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
  res <- .compute_estimable_contrast(m, linfct, confint = 0.95)
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
  res <- .compute_estimable_contrast(m, linfct, confint = 0.9)

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

test_that("compute_contrast estimates every row in the row space of a rank-deficient design", {
  withr::local_seed(7)
  d <- expand.grid(rep = 1:3, Treatment = c("A", "B"), Background = c("X", "Z"))
  d <- d[!(d$Treatment == "A" & d$Background == "X"), ]
  cell <- paste0(d$Treatment, d$Background)
  d$y <- stats::rnorm(nrow(d), mean = c(AZ = 10, BX = 12, BZ = 15)[cell])
  m <- lm(y ~ Treatment * Background, data = d)
  expect_true(is.na(coefficients(m)[["TreatmentB:BackgroundZ"]]))
  # coefficient weights of cell differences in the complete design: BZ - BX, BZ - AZ and AZ - AX (AX empty)
  linfct <- rbind(background_in_B = c(0, 0, 1, 1), treatment_in_Z = c(0, 1, 0, 1), background_in_A = c(0, 0, 1, 0))
  colnames(linfct) <- names(coefficients(m))
  res <- .compute_estimable_contrast(m, linfct)

  means <- tapply(d$y, cell, mean)
  expect_equal(res$estimate[1:2], unname(c(means["BZ"] - means["BX"], means["BZ"] - means["AZ"])))
  # the model is saturated on the observed cells: a difference of two cell means of 3 has SE sigma * sqrt(2 / 3)
  expect_equal(res$std.error[1:2], rep(sigma(m) * sqrt(2 / 3), 2))
  expect_true(all(is.na(res[3, c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")])))
})

test_that("compute_contrast follows the weighted design of an rlm fit whose IRLS weights alias a coefficient", {
  # one protein of the IonStar benchmark: the robust fit weights the two `a` observations near 0
  d <- data.frame(
    y = c(22.41270, 22.36527, 21.48429, 22.19423, 22.64496),
    dilution = c("b", "c", "d", "a", "a")
  )
  m <- MASS::rlm(y ~ dilution, data = d, method = "M")
  expect_true(is.na(stats::coef(m)[["dilutiond"]]))
  expect_false("dilutiond" %in% colnames(stats::vcov(m)))
  linfct <- rbind(c_vs_b = c(0, -1, 1, 0), d_vs_c = c(0, 0, -1, 1), b_vs_a = c(0, 1, 0, 0))
  colnames(linfct) <- names(stats::coef(m))
  res <- .compute_estimable_contrast(m, linfct)

  # b, c and d keep their full weight: their differences are those of the single observations
  expect_equal(res$estimate[1:2], c(22.36527 - 22.41270, 21.48429 - 22.36527), tolerance = 1e-6)
  # the `a` observations carry no weight, so contrasts with `a` are not estimable
  expect_true(is.na(res$estimate[3]))
})

test_that("linfct_matrix_contrasts names an unnamed contrast contrast_<i>", {
  linfct <- .linfct_from_model(.sim_make_model_lm("parallel3"), as_list = FALSE)
  res <- expect_no_warning(.linfct_matrix_contrasts(linfct, c("group_A - group_Ctrl", AB = "group_A - group_B")))
  expect_identical(rownames(res), c("contrast_1", "AB"))
  expect_equal(res["contrast_1", ], linfct["group_A", ] - linfct["group_Ctrl", ])
})

test_that("an infinite variance of an unused coefficient does not reach the contrast (rfit on small groups)", {
  linfct <- rbind(uses_b = c(0, 1), uses_a = c(1, 0))
  V <- matrix(c(Inf, 0, 0, 2), 2)
  expect_equal(prolfqua:::.contrast_variance(linfct, V), c(2, Inf))
})

test_that("model_coeff_matrix matches factor levels with non-word characters and 3-factor interactions exactly", {
  # Level "A" must not match "A.1" or "A-1" via regex word boundaries, and "+" must be matched literally
  df <- data.frame(
    protein_Id = rep("P1", 12),
    sample = paste0("S", 1:12),
    group_ = factor(rep(c("Ctrl", "A", "A.1", "A+B"), each = 3), levels = c("Ctrl", "A", "A.1", "A+B")),
    intensity = c(10, 10.1, 9.9, 15, 15.1, 14.9, 12, 12.1, 11.9, 18, 18.1, 17.9),
    nr_children = 3,
    identification_confidence = 1
  )
  cfg <- AnalysisConfiguration$new()
  cfg$file_name <- "sample"
  cfg$hierarchy <- list(protein_Id = "protein_Id")
  cfg$factors <- list(group_ = "group_")
  cfg$work_intensity <- "intensity"
  cfg$nr_children <- "nr_children"
  cfg$identification_confidence <- "identification_confidence"
  cfg$is_response_transformed <- TRUE
  lfq <- LFQData$new(setup_analysis(df, cfg), cfg)

  contrasts <- c(
    A_vs_Ctrl = "group_A - group_Ctrl",
    A1_vs_Ctrl = "group_A.1 - group_Ctrl",
    AB_vs_Ctrl = "`group_A+B` - group_Ctrl"
  )
  res_lm <- build_contrast_analysis(lfq, "~ group_", contrasts, method = "lm")$get_contrasts()
  expect_equal(res_lm$diff[match(names(contrasts), res_lm$contrast)], c(5, 2, 8), tolerance = 1e-6)

  res_limma <- build_contrast_analysis(lfq, "~ group_", contrasts, method = "limma")$get_contrasts()
  expect_equal(res_limma$diff[match(names(contrasts), res_limma$contrast)], c(5, 2, 8), tolerance = 1e-6)

  # 3-factor model: non-adjacent interaction term Aa1:Cc1 must apply to Aa1:Bb0:Cc1 and Aa1:Bb1:Cc1
  withr::local_seed(1)
  d3 <- expand.grid(rep = 1:3, A = c("a0", "a1"), B = c("b0", "b1"), C = c("c0", "c1"))
  d3$y <- with(
    d3,
    10 + 2 * (A == "a1") + 3 * (B == "b1") + 4 * (C == "c1") + 5 * (A == "a1" & C == "c1") + rnorm(nrow(d3), 0, 0.01)
  )
  m3 <- lm(y ~ A * B * C, data = d3)
  cm <- .model_coeff_matrix(m3)$coeff_matrix
  expected_ac <- as.numeric(rownames(cm) %in% c("Aa1:Bb0:Cc1", "Aa1:Bb1:Cc1"))
  expect_equal(unname(cm[, "Aa1:Cc1"]), expected_ac)
})
