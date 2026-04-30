test_that("fitFDistRobustly_LG preserves non-covariate robust fit behavior", {
  set.seed(1)
  x <- exp(rnorm(20, 0, 0.5))
  df1 <- sample(3:8, 20, replace = TRUE)

  fit <- prolfqua:::fitFDistRobustly_LG(x, df1)

  expect_equal(fit$scale, 1.202628, tolerance = 1e-6)
  expect_equal(fit$df2, Inf)
  expect_equal(fit$df2.shrunk, rep(Inf, 20))
})

test_that("fitFDistRobustly_LG preserves covariate trend behavior", {
  set.seed(1)
  x <- exp(rnorm(20, 0, 0.5))
  df1 <- sample(3:8, 20, replace = TRUE)

  fit <- prolfqua:::fitFDistRobustly_LG(x, df1, covariate = seq_along(x))

  expect_equal(
    fit$scale,
    c(
      0.9457446,
      1.0800594,
      1.2248818,
      1.3698919,
      1.5003645,
      1.5980796,
      1.6438567,
      1.6233316,
      1.5479917,
      1.4447611,
      1.3378039,
      1.2458362,
      1.1827813,
      1.1603041,
      1.1832065,
      1.2471774,
      1.3498923,
      1.4903815,
      1.6674290,
      1.8779027
    ),
    tolerance = 1e-6
  )
  expect_equal(fit$df2, Inf)
  expect_equal(fit$df2.shrunk, rep(Inf, 20))
})

test_that("fitFDistRobustly_LG handles very small variances", {
  set.seed(1)
  x <- exp(rnorm(20, 0, 0.5))
  x[1] <- 1e-20
  df1 <- sample(3:8, 20, replace = TRUE)

  expect_warning(
    fit <- prolfqua:::fitFDistRobustly_LG(x, df1),
    "One very small variance detected"
  )

  expect_equal(fit$scale, 1.314054, tolerance = 1e-6)
  expect_equal(fit$df2, 9.88225, tolerance = 1e-5)
  expect_equal(fit$df2.outlier, 2.50479, tolerance = 1e-5)
  expect_true(all(is.finite(fit$df2.shrunk)))
  expect_equal(length(fit$tail.p.value), length(x))
  expect_equal(length(fit$prob.outlier), length(x))
})

test_that("fitFDistRobustly_LG validates lengths", {
  expect_error(
    prolfqua:::fitFDistRobustly_LG(1:3, df1 = 1:2),
    "x and df1 are different lengths"
  )
  expect_error(
    prolfqua:::fitFDistRobustly_LG(1:3, df1 = 2, covariate = 1:2),
    "x and covariate are different lengths"
  )
})

test_that("fitFDistRobustly_LG reconstructs results after missing inputs", {
  fit <- prolfqua:::fitFDistRobustly_LG(
    x = c(1, NA, 2, 3, 4),
    df1 = c(3, 0, 4, 5, 6)
  )

  expect_equal(fit$scale, 2.5)
  expect_equal(fit$df2, Inf)
  expect_equal(fit$df2.shrunk, rep(Inf, 5))
})

test_that("squeezeVarRob uses robust prior degrees of freedom", {
  set.seed(1)
  x <- exp(rnorm(20, 0, 0.5))
  df1 <- sample(3:8, 20, replace = TRUE)

  fit <- prolfqua::squeezeVarRob(x, df1, robust = TRUE)

  expect_equal(fit$df.prior, rep(Inf, 20))
  expect_equal(fit$var.prior, 1.202628, tolerance = 1e-6)
  expect_equal(fit$var.post, rep(1.202628, 20), tolerance = 1e-6)
})
