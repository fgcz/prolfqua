test_that("moderated_p_limma moderates via limma::squeezeVar", {
  set.seed(42)
  n <- 60
  sigma <- sqrt(rexp(n))
  diff <- rnorm(n)
  std_error <- sigma * runif(n, 0.4, 1.4)
  contrast_df <- data.frame(
    contrast = "A_vs_B",
    sigma = sigma,
    df = sample(3:9, n, replace = TRUE), # unequal df exercises the modern estimator
    statistic = diff / std_error,
    std.error = std_error,
    diff = diff
  )
  res <- prolfqua::moderated_p_limma(contrast_df)
  sv <- limma::squeezeVar(contrast_df$sigma^2, df = contrast_df$df)

  expect_equal(res$moderated.var.post, sv$var.post, tolerance = 1e-10)
  expect_equal(unique(res$moderated.df.prior), sv$df.prior, tolerance = 1e-10)
})

test_that("moderated_p_limma honors custom columns and confidence levels", {
  contrast_df <- data.frame(
    sigma = c(0.25, 0.32, 0.28, 0.40, 0.35),
    residual_df = rep(6, 5),
    statistic = c(2.1, -1.8, 0.5, 3.0, -2.2),
    std.error = c(0.8 / 2.1, 0.6 / 1.8, 0.2 / 0.5, 1.2 / 3.0, 0.9 / 2.2),
    effect = c(0.8, -0.6, 0.2, 1.2, -0.9)
  )

  result <- moderated_p_limma(
    contrast_df,
    df = "residual_df",
    estimate = "effect",
    confint = 0.8
  )
  quantile <- -stats::qt(
    0.1,
    df = result$moderated.df.total
  )

  expect_equal(
    result$moderated.df.total,
    contrast_df$residual_df + result$moderated.df.prior
  )
  expect_equal(
    result$moderated.std.error,
    contrast_df$std.error * sqrt(result$moderated.var.post) / contrast_df$sigma
  )
  expect_equal(
    result$moderated.conf.low,
    contrast_df$effect -
      quantile * result$moderated.std.error
  )
  expect_equal(
    result$moderated.conf.high,
    contrast_df$effect +
      quantile * result$moderated.std.error
  )
})

test_that("moderated_p_limma rejects invalid variance floors", {
  contrast_df <- data.frame(
    sigma = c(0.25, 0.32, 0.28, 0.40, 0.35),
    df = rep(6, 5),
    statistic = c(2.1, -1.8, 0.5, 3.0, -2.2),
    std.error = c(0.8 / 2.1, 0.6 / 1.8, 0.2 / 0.5, 1.2 / 3.0, 0.9 / 2.2),
    diff = c(0.8, -0.6, 0.2, 1.2, -0.9)
  )
  invalid_floors <- list(0, -1, c(1, 2), NA_real_, Inf, "1")

  purrr::walk(
    invalid_floors,
    ~ expect_error(
      moderated_p_limma(
        contrast_df,
        variance_floor = .x
      ),
      "positive number",
      fixed = TRUE
    )
  )
})
