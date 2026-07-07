# Test harness for Items 5 & 6 vectorization
#
# Strategy: run original and _vectorized side by side, verify identical output.
# The originals are untouched and known to work; the _vectorized variants must
# produce the same results.

# ---- Section A: linfct_matrix_contrasts vs linfct_matrix_contrasts_vectorized ----

test_that("linfct_matrix_contrasts_vectorized matches original: parallel3", {
  set.seed(42)
  m <- sim_make_model_lm("parallel3")
  linfct <- linfct_from_model(m, as_list = FALSE)

  Contr <- c(
    "AvsCtrl" = "group_A - group_Ctrl",
    "BvsCtrl" = "group_B - group_Ctrl"
  )

  orig <- linfct_matrix_contrasts(linfct, Contr)
  vec <- linfct_matrix_contrasts_vectorized(linfct, Contr)

  expect_identical(dim(orig), dim(vec))
  expect_identical(rownames(orig), rownames(vec))
  expect_identical(colnames(orig), colnames(vec))
  expect_equal(orig, vec, tolerance = 1e-12)
})

test_that("linfct_matrix_contrasts_vectorized matches original: self-referencing (factors)", {
  set.seed(42)
  m <- sim_make_model_lm("factors")
  linfct <- linfct_from_model(m, as_list = FALSE)

  Contr <- c(
    "TreatmentA_vs_B" = "TreatmentA - TreatmentB",
    "BackgroundX_vs_Z" = "BackgroundX - BackgroundZ",
    "IntoflintoA" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
    "IntoflintoB" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`",
    "IntoflintoX" = "`TreatmentA:BackgroundX` - `TreatmentB:BackgroundX`",
    "IntoflintoZ" = "`TreatmentA:BackgroundZ` - `TreatmentB:BackgroundZ`",
    "interactXZ" = "IntoflintoX - IntoflintoZ",
    "interactAB" = "IntoflintoA - IntoflintoB"
  )

  orig <- linfct_matrix_contrasts(linfct, Contr)
  vec <- linfct_matrix_contrasts_vectorized(linfct, Contr)

  expect_equal(orig, vec, tolerance = 1e-12)

  # Sanity: interaction contrasts sum to 0 for additive model
  expect_equal(sum(vec["interactXZ", ]), 0)
  expect_equal(sum(vec["interactAB", ]), 0)
})

test_that("linfct_matrix_contrasts_vectorized matches original: interaction model", {
  set.seed(42)
  m <- sim_make_model_lm("interaction")
  linfct <- linfct_from_model(m, as_list = FALSE)

  Contr <- c(
    "TreatmentA_vs_B" = "TreatmentA - TreatmentB",
    "BackgroundX_vs_Z" = "BackgroundX - BackgroundZ",
    "IntoflintoA" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
    "IntoflintoB" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`",
    "IntoflintoX" = "`TreatmentA:BackgroundX` - `TreatmentB:BackgroundX`",
    "IntoflintoZ" = "`TreatmentA:BackgroundZ` - `TreatmentB:BackgroundZ`",
    "interactXZ" = "IntoflintoX - IntoflintoZ",
    "interactAB" = "IntoflintoA - IntoflintoB"
  )

  orig <- linfct_matrix_contrasts(linfct, Contr)
  vec <- linfct_matrix_contrasts_vectorized(linfct, Contr)

  expect_equal(orig, vec, tolerance = 1e-12)
  expect_equal(sum(vec["interactXZ", ]), 1)
  expect_equal(sum(vec["interactAB", ]), 1)
})

test_that("linfct_matrix_contrasts_vectorized matches original: single contrast", {
  set.seed(42)
  m <- sim_make_model_lm("parallel3")
  linfct <- linfct_from_model(m, as_list = FALSE)

  Contr <- c("AvsCtrl" = "group_A - group_Ctrl")

  orig <- linfct_matrix_contrasts(linfct, Contr)
  vec <- linfct_matrix_contrasts_vectorized(linfct, Contr)

  expect_equal(orig, vec, tolerance = 1e-12)
})

test_that("linfct_matrix_contrasts_vectorized matches original: invalid contrast warning", {
  set.seed(42)
  m <- sim_make_model_lm("parallel3")
  linfct <- linfct_from_model(m, as_list = FALSE)

  Contr <- c(
    "valid" = "group_A - group_Ctrl",
    "invalid" = "nonexistent_group - group_Ctrl"
  )

  expect_warning(
    orig <- linfct_matrix_contrasts(linfct, Contr),
    "linfct_matrix_contrasts"
  )
  expect_warning(
    vec <- linfct_matrix_contrasts_vectorized(linfct, Contr),
    "linfct_matrix_contrasts"
  )

  expect_equal(orig, vec, tolerance = 1e-12)
})

# ---- Section B: compute_contrast vs compute_contrast_vectorized ----

test_that("compute_contrast_vectorized matches original: full factors model", {
  set.seed(42)
  m <- sim_make_model_lm("factors")
  linfct <- linfct_from_model(m)$linfct_factors

  orig <- compute_contrast(m, linfct, confint = 0.95)
  vec <- compute_contrast_vectorized(m, linfct, confint = 0.95)

  expect_equal(nrow(orig), nrow(vec))
  expect_equal(orig$lhs, vec$lhs)
  expect_equal(orig$estimate, vec$estimate, tolerance = 1e-12)
  expect_equal(orig$std.error, vec$std.error, tolerance = 1e-12)
  expect_equal(orig$statistic, vec$statistic, tolerance = 1e-12)
  expect_equal(orig$p.value, vec$p.value, tolerance = 1e-12)
  expect_equal(orig$conf.low, vec$conf.low, tolerance = 1e-12)
  expect_equal(orig$conf.high, vec$conf.high, tolerance = 1e-12)
})

test_that("compute_contrast_vectorized matches original: parallel3 model", {
  set.seed(42)
  m <- sim_make_model_lm("parallel3")
  linfct <- linfct_from_model(m)$linfct_factors

  orig <- compute_contrast(m, linfct, confint = 0.95)
  vec <- compute_contrast_vectorized(m, linfct, confint = 0.95)

  expect_equal(orig$estimate, vec$estimate, tolerance = 1e-12)
  expect_equal(orig$std.error, vec$std.error, tolerance = 1e-12)
  expect_equal(orig$p.value, vec$p.value, tolerance = 1e-12)
})

test_that("compute_contrast_vectorized matches original: single-row linfct", {
  set.seed(42)
  m <- sim_make_model_lm("parallel3")
  linfct <- linfct_from_model(m)$linfct_factors
  linfct_one <- linfct[1, , drop = FALSE]

  orig <- compute_contrast(m, linfct_one, confint = 0.95)
  vec <- compute_contrast_vectorized(m, linfct_one, confint = 0.95)

  expect_equal(orig$estimate, vec$estimate, tolerance = 1e-12)
  expect_equal(orig$p.value, vec$p.value, tolerance = 1e-12)
  expect_equal(orig$lhs, vec$lhs)
})

test_that("compute_contrast_vectorized matches original: all pairwise contrasts", {
  set.seed(42)
  m <- sim_make_model_lm("parallel3")
  linfct_list <- linfct_from_model(m)
  all_pairs <- linfct_all_possible_contrasts(linfct_list$linfct_factors)

  orig <- compute_contrast(m, all_pairs, confint = 0.95)
  vec <- compute_contrast_vectorized(m, all_pairs, confint = 0.95)

  expect_equal(orig$lhs, vec$lhs)
  expect_equal(orig$estimate, vec$estimate, tolerance = 1e-12)
  expect_equal(orig$std.error, vec$std.error, tolerance = 1e-12)
  expect_equal(orig$p.value, vec$p.value, tolerance = 1e-12)
})

test_that("compute_contrast_vectorized matches original: different confint", {
  set.seed(42)
  m <- sim_make_model_lm("factors")
  linfct <- linfct_from_model(m)$linfct_factors

  orig_95 <- compute_contrast(m, linfct, confint = 0.95)
  vec_95 <- compute_contrast_vectorized(m, linfct, confint = 0.95)
  orig_99 <- compute_contrast(m, linfct, confint = 0.99)
  vec_99 <- compute_contrast_vectorized(m, linfct, confint = 0.99)

  expect_equal(orig_95$conf.low, vec_95$conf.low, tolerance = 1e-12)
  expect_equal(orig_95$conf.high, vec_95$conf.high, tolerance = 1e-12)
  expect_equal(orig_99$conf.low, vec_99$conf.low, tolerance = 1e-12)
  expect_equal(orig_99$conf.high, vec_99$conf.high, tolerance = 1e-12)
})

test_that("compute_contrast_vectorized matches original: incomplete model (NA coefficients)", {
  set.seed(42)
  df <- data.frame(
    y = c(rnorm(3, mean = 10), rnorm(3, mean = 12)),
    group = factor(c("A", "A", "A", "B", "B", "B"), levels = c("A", "B", "C"))
  )
  m <- lm(y ~ group, data = df)
  expect_true(is.na(coefficients(m)["groupC"]))

  linfct <- matrix(0, nrow = 2, ncol = 3, dimnames = list(c("AvsB", "AvsC"), c("(Intercept)", "groupB", "groupC")))
  linfct["AvsB", ] <- c(0, -1, 0)
  linfct["AvsC", ] <- c(0, 0, -1)

  orig <- compute_contrast(m, linfct, confint = 0.95)
  vec <- compute_contrast_vectorized(m, linfct, confint = 0.95)

  # Both should have same NA pattern
  expect_equal(is.na(orig$estimate), is.na(vec$estimate))
  expect_equal(is.na(orig$std.error), is.na(vec$std.error))
  expect_equal(is.na(orig$p.value), is.na(vec$p.value))

  # Valid rows must match numerically
  expect_equal(orig$estimate[1], vec$estimate[1], tolerance = 1e-12)
  expect_equal(orig$std.error[1], vec$std.error[1], tolerance = 1e-12)
  expect_equal(orig$p.value[1], vec$p.value[1], tolerance = 1e-12)

  # Row order preserved
  expect_equal(orig$lhs, vec$lhs)
})

test_that("compute_contrast_vectorized matches original: simulated incomplete models", {
  set.seed(42)
  mod <- sim_build_models_lm("parallel3", Nprot = 50, with_missing = TRUE, weight_missing = 2)

  has_na <- vapply(
    mod$model_df$linear_model,
    function(m) {
      inherits(m, "lm") && any(is.na(coefficients(m)))
    },
    logical(1)
  )

  skip_if(!any(has_na), "No models with NA coefficients in simulated data")

  idx <- which(has_na)[1]
  m <- mod$model_df$linear_model[[idx]]

  linfct_list <- linfct_from_model(m)
  all_pairs <- linfct_all_possible_contrasts(linfct_list$linfct_factors)

  orig <- compute_contrast(m, all_pairs, confint = 0.95)
  vec <- compute_contrast_vectorized(m, all_pairs, confint = 0.95)

  # Same NA pattern
  expect_equal(is.na(orig$estimate), is.na(vec$estimate))

  # Valid rows match numerically
  valid <- !is.na(orig$estimate)
  if (any(valid)) {
    expect_equal(orig$estimate[valid], vec$estimate[valid], tolerance = 1e-12)
    expect_equal(orig$std.error[valid], vec$std.error[valid], tolerance = 1e-12)
    expect_equal(orig$p.value[valid], vec$p.value[valid], tolerance = 1e-12)
  }

  expect_equal(orig$lhs, vec$lhs)
})

test_that("compute_contrast_vectorized matches original: rank-deficient model, canceling NA weights", {
  # Deterministic rank-deficient fit: B is perfectly aliased with A, so coefficients
  # Bb2/Bb3 are NA. A contrast placing +1/-1 on those two NA coefficients has a signed
  # weight sum of 0 -- the previous vectorized guard wrongly treated it as valid and
  # returned estimate = 0 / p = NaN, while the loop path correctly returns NA.
  d <- data.frame(
    y = c(1.2, -0.3, 0.8, 0.1, 2.1, -1.0, 0.4, 0.9, -0.6, 1.5, 0.2, -0.4),
    A = factor(rep(c("a1", "a2", "a3"), each = 4)),
    B = factor(rep(c("b1", "b2", "b3"), each = 4))
  )
  m <- lm(y ~ A + B, data = d)
  expect_true(any(is.na(coefficients(m)))) # confirm rank deficiency

  na_names <- names(which(is.na(coefficients(m))))
  expect_gte(length(na_names), 2)
  cn <- names(coefficients(m))
  linfct <- matrix(0, nrow = 1, ncol = length(cn), dimnames = list("cancel", cn))
  linfct[1, na_names[1]] <- 1
  linfct[1, na_names[2]] <- -1

  orig <- compute_contrast(m, linfct, confint = 0.95)
  vec <- compute_contrast_vectorized(m, linfct, confint = 0.95)

  # Loop path refuses the contrast (references non-estimable coefficients).
  expect_true(is.na(orig$estimate))
  # Vectorized path must agree -- not a spurious 0 estimate.
  expect_equal(is.na(orig$estimate), is.na(vec$estimate))
  expect_equal(is.na(orig$p.value), is.na(vec$p.value))
  expect_true(is.na(vec$estimate))
})

# ---- Section C: options(prolfqua.vectorize) hot-swap through facade ----

test_that("prolfqua.vectorize option dispatches correctly through build_contrast_analysis", {
  istar <- sim_lfq_data_protein_config(Nprot = 20)
  lfqdata <- LFQData$new(istar$data, istar$config)
  lfqdata <- lfqdata$get_Transformer()$log2()$lfq
  contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")

  # Original path (default)
  options(prolfqua.vectorize = FALSE)
  fa_orig <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")
  res_orig <- fa_orig$get_contrasts()

  # Vectorized path
  options(prolfqua.vectorize = TRUE)
  fa_vec <- build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm")
  res_vec <- fa_vec$get_contrasts()

  # Restore default
  options(prolfqua.vectorize = FALSE)

  expect_equal(nrow(res_orig), nrow(res_vec))
  expect_equal(res_orig$diff, res_vec$diff, tolerance = 1e-12)
  expect_equal(res_orig$std.error, res_vec$std.error, tolerance = 1e-12)
  expect_equal(res_orig$p.value, res_vec$p.value, tolerance = 1e-12)
  expect_equal(res_orig$FDR, res_vec$FDR, tolerance = 1e-12)
})
