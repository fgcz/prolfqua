test_that("moderated_p_limma moderates via limma::squeezeVar", {
  set.seed(42)
  n <- 60
  contrast_df <- data.frame(
    contrast = "A_vs_B",
    sigma = sqrt(rexp(n)),
    df = sample(3:9, n, replace = TRUE), # unequal df exercises the modern estimator
    statistic = rnorm(n),
    diff = rnorm(n)
  )
  res <- prolfqua::moderated_p_limma(contrast_df)
  sv <- limma::squeezeVar(contrast_df$sigma^2, df = contrast_df$df)

  expect_equal(res$moderated.var.post, sv$var.post, tolerance = 1e-10)
  expect_equal(unique(res$moderated.df.prior), sv$df.prior, tolerance = 1e-10)
})
