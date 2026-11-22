test_that("FDR planning reproduces the documented standardized-effect example", {
  alpha <- fdr_sig_level(fdr = 0.05, pi0 = 0.9, power = 0.8)
  expect_equal(alpha, 0.00467836257309942)
  expected <- stats::power.t.test(delta = 1, sd = 1, sig.level = alpha, power = 0.8)
  expect_equal(ceiling(expected$n), 29)
  expect_equal(.power_t_test_n(1, 1, 0.8, alpha, 2)$N, 29)
})

test_that("FDR probability inputs have classed errors", {
  invalid <- list(0, 1, -0.1, Inf, NA_real_, NaN, NULL, c(0.1, 0.2), "0.05", TRUE, 1i, matrix(0.1))
  for (argument in c("fdr", "pi0", "power")) {
    for (value in invalid) {
      expect_error(
        do.call(fdr_sig_level, setNames(list(value), argument)),
        argument,
        class = "prolfqua_error_bad_argument"
      )
    }
  }
  expect_error(fdr_sig_level(fdr = 0.9, pi0 = 0.01), "per-test significance", class = "prolfqua_error_bad_argument")
})

test_that("pilot prior inputs are validated before fitting", {
  valid <- list(var = rep(1, 20), df = 8, delta_mean = 1, max_n = 10)
  invalid <- list(
    var = list("1", list(1), matrix(1, 20), c(1, 2, -1), complex(20), rep(0, 20), c(0, 0, 1)),
    df = list("8", c(4, 8), matrix(8, 20), complex(20)),
    delta_mean = list(NA_real_, Inf, numeric(), c(0.1, 0.2), "1", 1i),
    delta_se = list(-0.1, Inf, NA_real_, c(0.1, 0.2), "1", 1i),
    max_n = list(1, 2.5, Inf, NA_real_, c(20, 30), "20", 1i)
  )
  for (argument in names(invalid)) {
    for (value in invalid[[argument]]) {
      args <- valid
      args[argument] <- list(value)
      expect_error(do.call(power_fdr_prior_n, args), argument, class = "prolfqua_error_bad_argument")
    }
  }
  for (variances in list(numeric(), 1, c(1, 2), rep(NA_real_, 4))) {
    expect_error(
      power_fdr_prior_n(variances, df = 8, delta_mean = 1),
      "at least three",
      class = "prolfqua_error_bad_argument"
    )
  }
  expect_error(power_fdr_prior_n(rep(1, 20), df = 0, delta_mean = 1), class = "prolfqua_error_bad_argument")
})

test_that("infinite prior df reproduces the fixed-variance t-test", {
  result <- power_fdr_prior_n(rep(1, 100), df = 8, delta_mean = 1, pi0 = 0.9, max_n = 60)
  expect_identical(result$df_prior, Inf)
  expect_equal(result$var_prior, 1)
  expect_identical(result$N, 29L)
  expected <- vapply(
    2:60,
    function(n) {
      stats::power.t.test(n = n, delta = 1, sd = 1, sig.level = result$sig_level, strict = TRUE)$power
    },
    numeric(1)
  )
  expect_equal(result$power_curve$power, expected, tolerance = 1e-12)
  expect_equal(result$power_achieved, expected[result$N - 1L])
  expect_gt(result$power_achieved, result$power)
  expect_lt(result$power_curve$power[result$N - 2L], result$power)
})

test_that("normal-effect quadrature agrees with analytic normal mixing at fixed variance", {
  n <- c(2, 10, 30, 60)
  alpha <- fdr_sig_level()
  scale <- sqrt(1 + 0.5^2 * n / (2 * 0.5))
  critical <- stats::qt(alpha / 2, df = 2 * n - 2, lower.tail = FALSE) / scale
  ncp <- sqrt(n / (2 * 0.5)) / scale
  expected <- stats::pt(critical, df = 2 * n - 2, ncp = ncp, lower.tail = FALSE) +
    stats::pt(-critical, df = 2 * n - 2, ncp = ncp)
  actual <- .fdr_prior_power(n, 0.5, Inf, 1, 0.5, alpha)
  expect_equal(actual, expected, tolerance = 2e-4)
})

test_that("unreachable targets retain achieved power at max_n and use both tails", {
  result <- power_fdr_prior_n(rep(1, 20), df = 8, delta_mean = 0, max_n = 10)
  expect_identical(result$N, NA_integer_)
  expect_equal(result$power_curve$n, 2:10)
  expect_equal(result$power_curve$power, rep(result$sig_level, 9))
  expect_equal(result$power_achieved, result$sig_level)
  positive <- power_fdr_prior_n(rep(1, 20), df = 8, delta_mean = 1, delta_se = 0.2, max_n = 10)
  negative <- power_fdr_prior_n(rep(1, 20), df = 8, delta_mean = -1, delta_se = 0.2, max_n = 10)
  expect_equal(positive$power_curve, negative$power_curve)
  small <- power_fdr_prior_n(rep(1, 20), df = 8, delta_mean = 1, max_n = 2)
  expect_equal(nrow(small$power_curve), 1L)
})

test_that("prior fitting filters unusable observations and preserves original inputs", {
  var <- c(0.1, 0.2, 0.4, 0.6, NA, Inf, NaN, 0.9, 0.5, 0.5, 0.3)
  df <- c(4, 6, 8, 10, 8, 8, 8, 0, NA, Inf, -1)
  actual <- power_fdr_prior_n(var, df, delta_mean = 1, max_n = 10)
  expected <- power_fdr_prior_n(var[1:4], df[1:4], delta_mean = 1, max_n = 10)
  expect_identical(actual$n_features_used, 4L)
  expect_equal(actual$var_prior, expected$var_prior)
  expect_equal(actual$df_prior, expected$df_prior)
  expect_equal(actual$power_curve, expected$power_curve)
  expect_identical(actual$var, var)
  expect_identical(actual$df, df)
})

test_that("pilot variances recover a known scaled inverse-chi-square prior", {
  withr::local_seed(42)
  var <- 0.25 * stats::rf(20000, df1 = 8, df2 = 12)
  result <- power_fdr_prior_n(var, df = 8, delta_mean = 0.8, max_n = 60)
  expect_lt(abs(result$var_prior / 0.25 - 1), 0.05)
  expect_lt(abs(result$df_prior / 12 - 1), 0.1)
  expect_identical(result$n_features_used, 20000L)
})

test_that("large finite prior df approaches fixed variance without unstable integration", {
  alpha <- fdr_sig_level()
  finite <- .fdr_prior_power(2:60, 0.5, 200, 1, 0.2, alpha)
  fixed <- .fdr_prior_power(2:60, 0.5, Inf, 1, 0.2, alpha)
  expect_true(all(is.finite(finite)))
  expect_lt(max(abs(finite - fixed)), 0.01)
  expect_true(all(diff(finite) >= 0))
})

test_that("small prior df quadrature agrees with integrated normal-mixture power", {
  alpha <- fdr_sig_level()
  n <- c(4, 20, 60)
  expected <- vapply(
    n,
    function(size) {
      stats::integrate(
        function(u) {
          precision <- stats::qgamma(u, shape = 1, rate = 0.5)
          scale <- sqrt(1 + 0.3^2 * size * precision / 2)
          critical <- stats::qt(alpha / 2, 2 * size - 2, lower.tail = FALSE) / scale
          ncp <- sqrt(size * precision / 2) / scale
          stats::pt(critical, 2 * size - 2, ncp, lower.tail = FALSE) + stats::pt(-critical, 2 * size - 2, ncp)
        },
        0,
        1,
        rel.tol = 1e-8
      )$value
    },
    numeric(1)
  )
  actual <- .fdr_prior_power(n, 0.5, 2, 1, 0.3, alpha)
  expect_equal(actual, expected, tolerance = 1e-4)
})

test_that("LFQData pilot planning pools observed within-group variances", {
  sim <- sim_lfq_data_protein_config(n_proteins = 20, seed = 42)
  lfq <- LFQData$new(sim$data, sim$config)$get_Transformer()$log2()$lfq
  pooled <- poolvar(summarize_stats(lfq), lfq$get_config())
  expected <- power_fdr_prior_n(pooled$var, pooled$df, delta_mean = 0.1, max_n = 30)
  actual <- lfq$get_Stats(stats = "all")$power_fdr_pilot(delta_mean = 0.1, max_n = 30)
  expect_equal(actual, expected)
  expect_equal(lfq$get_Stats(stats = "interaction")$power_fdr_pilot(0.1, max_n = 30), expected)
  expect_equal(sum(pooled$df), sum(pmax(summarize_stats(lfq)$nrMeasured - 1, 0)))
  data <- lfq$data_long()
  first <- !duplicated(data[, c(lfq$hierarchy_keys(), lfq$relevant_factor_keys())])
  lfq$set_data(data[first, ])
  expect_error(lfq$get_Stats()$power_fdr_pilot(0.1), "needs replicates", class = "prolfqua_error_bad_argument")
})

fdr_reference <- function(fun, ...) {
  grDevices::pdf(NULL)
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device), add = TRUE)
  muffle_warnings_matching(fun(...), "full precision may not have been achieved")
}

test_that("fixed-effect FDR sample sizes agree with ssize.fdr on a grid", {
  skip_if_not_installed("ssize.fdr")
  grid <- expand.grid(delta = c(0.5, 1, 2), pi0 = c(0.8, 0.9, 0.99), fdr = c(0.05, 0.1), power = c(0.8, 0.9))
  for (i in seq_len(nrow(grid))) {
    x <- grid[i, ]
    alpha <- fdr_sig_level(x$fdr, x$pi0, x$power)
    n <- ceiling(stats::power.t.test(delta = x$delta, sd = 1, sig.level = alpha, power = x$power)$n)
    reference <- fdr_reference(
      ssize.fdr::ssize.twoSamp,
      delta = x$delta,
      sigma = 1,
      fdr = x$fdr,
      pi0 = x$pi0,
      power = x$power,
      maxN = n + 1
    )
    expect_equal(n, unname(reference$ssize[1, "ssize"]))
  }
})

test_that("variance-prior FDR sample sizes agree with stable ssize.fdr cases", {
  skip_if_not_installed("ssize.fdr")
  grid <- expand.grid(df_prior = c(4, 12, 40), delta_se = c(0, 0.3), pi0 = c(0.9, 0.95))
  for (i in seq_len(nrow(grid))) {
    x <- grid[i, ]
    var <- 0.5 * stats::qf((seq_len(5000) - 0.5) / 5000, df1 = 8, df2 = x$df_prior)
    actual <- power_fdr_prior_n(var, df = 8, delta_mean = 1, delta_se = x$delta_se, pi0 = x$pi0, max_n = 60)
    reference <- fdr_reference(
      ssize.fdr::ssize.twoSampVary,
      deltaMean = 1,
      deltaSE = x$delta_se,
      a = actual$df_prior / 2,
      b = actual$df_prior * actual$var_prior / 2,
      pi0 = x$pi0,
      maxN = 60
    )
    expect_lte(abs(actual$N - unname(reference$ssize[1, "ssize"])), 1)
  }
})

test_that("power methods error on untransformed data", {
  x <- sim_lfq_data_protein_config(n_proteins = 5, seed = 42)
  stats <- LFQData$new(x$data, x$config)$get_Stats()
  expect_error(stats$power_t_test_quantiles(), "power_t_test_quantiles", class = "prolfqua_error_bad_argument")
  expect_error(stats$power_t_test(), "power_t_test", class = "prolfqua_error_bad_argument")
  expect_error(stats$power_fdr_pilot(delta_mean = 0.1), "power_fdr_pilot", class = "prolfqua_error_bad_argument")
})

test_that("lfq_power_t_test_proteins gives the power.t.test sample size, floored at min.n", {
  stats_res <- data.frame(protein_Id = c("tiny", "unit"), sd = c(0.01, 1))
  res <- lfq_power_t_test_proteins(stats_res, delta = 1, power = 0.8, sig.level = 0.05, min.n = 1.5)
  # n for delta = sd = 1 at 80% power and 5% two-sided level
  expect_equal(res$N_exact[res$protein_Id == "unit"], 16.71476814, tolerance = 1e-8)
  expect_equal(res$N[res$protein_Id == "unit"], 17)
  # a standard deviation below the floor gives the floor sample size
  expect_equal(res$N_exact[res$protein_Id == "tiny"], 1.5, tolerance = 1e-4)
  expect_equal(res$N[res$protein_Id == "tiny"], 2)
  # one row per protein and delta
  expect_equal(nrow(lfq_power_t_test_proteins(stats_res, delta = c(0.5, 1, 2))), 6)
})

test_that("LFQDataStats power methods equal power.t.test on the protein and quantile standard deviations", {
  x <- sim_lfq_data_protein_config(n_proteins = 20, with_missing = FALSE, seed = 42)
  lfq <- LFQData$new(x$data, x$config)
  lfq <- lfq$get_Transformer()$log2()$lfq
  stats <- lfq$get_Stats()

  per_protein <- stats$power_t_test(delta = 1)
  row <- per_protein[which.max(per_protein$sd), ]
  expect_equal(row$N_exact, stats::power.t.test(delta = 1, sd = row$sd, power = 0.8)$n, tolerance = 1e-6)
  expect_equal(nrow(per_protein), nrow(stats::na.omit(stats$stats())))

  quantiles <- stats$power_t_test_quantiles(probs = c(0.5, 0.9), delta = 1)
  expect_true(all(quantiles$probs %in% c(0.5, 0.9)))
  q <- quantiles[which.max(quantiles$quantiles), ]
  expect_equal(
    unname(q$N_exact),
    stats::power.t.test(delta = 1, sd = unname(q$quantiles), power = 0.8)$n,
    tolerance = 1e-6
  )
  # larger standard deviations need more samples
  expect_true(all(diff(quantiles$N_exact[order(quantiles$quantiles)]) >= 0))
})
