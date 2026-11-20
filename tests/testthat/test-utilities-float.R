# Floating-point comparisons use tolerances (srr G3.0).

test_that(".is_zero uses a tolerance and treats NA as non-zero", {
  expect_identical(.is_zero(c(0, 1e-12, -1e-12, 1e-6, NA)), c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("a near-zero ROPECA statistic gives NA standard error, not a huge one", {
  res <- data.frame(
    beta.based.significance = c(0.5, 0.01),
    FDR.beta.based.significance = c(0.5, 0.02),
    statistic = c(1e-14, 2),
    diff = c(0.3, 1),
    mad.estimate = 0.2,
    n_not_na = 3,
    protein_Id = c("a", "b"),
    modelName = "ROPECA",
    estimate_type = "observed",
    contrast = "AvsC",
    avgAbd = 10
  )
  out <- prolfqua:::.ropeca_standard_columns(res, "protein_Id")
  expect_true(is.na(out$std.error[1]))
  expect_equal(out$std.error[2], 0.5)
})

test_that("near-zero sample MADs are treated as zero when scaling", {
  data <- matrix(c(1, 1 + 1e-12, 1, 2, 3, 4), nrow = 3)
  scales <- list(medians = c(1, 3), mads = c(1e-13, 1))
  expect_warning(
    out <- prolfqua:::.apply_robscales(data, scales, preserve_mean = FALSE, dim = 2, zero_mad_warning = "zero MAD"),
    "zero MAD"
  )
  expect_equal(out[, 2], c(-1, 0, 1))
})
