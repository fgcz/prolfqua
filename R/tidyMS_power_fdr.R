#' Per-test significance level for FDR-controlled sample-size planning
#'
#' Convert a target false discovery rate and average power into the per-test
#' significance level used by an ordinary two-sided t-test. Pass the result to
#' [stats::power.t.test()], [lfq_power_t_test_proteins()], or
#' [lfq_power_t_test_quantiles()].
#'
#' @param fdr one target false discovery rate, strictly between zero and one.
#' @param pi0 one assumed proportion of truly null features, strictly between
#'   zero and one; it is not estimated from the pilot.
#' @param power one target average power among non-null features, strictly
#'   between zero and one.
#' @return The scalar per-test significance level
#'   \eqn{\alpha_* = power \times fdr \times (1-pi0) / (pi0 \times (1-fdr))}.
#'   Combinations giving a significance level outside (0, 1) are rejected.
#' @details The Jung/Liu--Hwang planning relationship assumes independent tests
#'   or an adequate approximation under weak dependence, and many tested
#'   features. It controls the expected false-discovery proportion under the
#'   assumed null fraction and effect distribution, not the realized FDR of
#'   every experiment. All sample sizes are per group in a balanced two-group
#'   design with equal within-group variances.
#' @references Jung SH (2005). Sample size for FDR-control in microarray data
#'   analysis. Bioinformatics 21:3097--3104. \doi{10.1093/bioinformatics/bti456}.
#'   Liu P, Hwang JTG (2007). Quick calculation for sample size while controlling
#'   false discovery rate with application to microarray analysis.
#'   Bioinformatics 23:739--746. \doi{10.1093/bioinformatics/btl664}.
#' @family stats
#' @export
#' @examples
#' alpha <- fdr_sig_level(fdr = 0.05, pi0 = 0.9, power = 0.8)
#' alpha
#' ceiling(stats::power.t.test(delta = 1, sd = 1, sig.level = alpha, power = 0.8)$n)
fdr_sig_level <- function(fdr = 0.05, pi0 = 0.95, power = 0.8) {
  .check_fdr_probability(fdr, "fdr")
  .check_fdr_probability(pi0, "pi0")
  .check_fdr_probability(power, "power")
  alpha <- power * fdr * (1 - pi0) / (pi0 * (1 - fdr))
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    abort_bad_argument("fdr, pi0, power", "give a per-test significance level strictly between zero and one")
  }
  alpha
}

#' Recognize a real numeric vector without matrix dimensions
#' @param x input value
#' @noRd
.fdr_numeric_vector <- function(x) is.numeric(x) && !is.complex(x) && is.null(dim(x))

#' Recognize a finite real numeric scalar
#' @param x input value
#' @noRd
.fdr_numeric_scalar <- function(x) .fdr_numeric_vector(x) && length(x) == 1L && is.finite(x)

#' Validate one FDR-planning probability
#' @param value probability value
#' @param argument argument name for the error
#' @noRd
.check_fdr_probability <- function(value, argument) {
  if (!.fdr_numeric_scalar(value) || value <= 0 || value >= 1) {
    abort_bad_argument(argument, "be one finite number strictly between zero and one")
  }
}

#' Validate effect-distribution and sample-size inputs
#' @param delta_mean planned mean response difference
#' @param delta_se planned standard deviation of differences among non-null features
#' @param max_n largest sample count per group
#' @noRd
.check_fdr_effect_inputs <- function(delta_mean, delta_se, max_n) {
  if (!.fdr_numeric_scalar(delta_mean)) {
    abort_bad_argument("delta_mean", "be one finite number")
  }
  if (!.fdr_numeric_scalar(delta_se) || delta_se < 0) {
    abort_bad_argument("delta_se", "be one finite non-negative number")
  }
  if (!.fdr_numeric_scalar(max_n) || max_n < 2 || max_n != floor(max_n)) {
    abort_bad_argument("max_n", "be one whole number of at least two samples per group")
  }
}

#' Select variance estimates with usable residual degrees of freedom
#' @param var per-feature variance estimates
#' @param df one residual df or a vector parallel to var
#' @noRd
.fdr_prior_inputs <- function(var, df) {
  if (!.fdr_numeric_vector(var)) {
    abort_bad_argument("var", "be a real numeric vector of variance estimates")
  }
  if (!.fdr_numeric_vector(df) || !length(df) %in% c(1L, length(var))) {
    abort_bad_argument("df", "be a real numeric vector of length one or length(var)")
  }
  if (any(is.finite(var) & var < 0)) {
    abort_bad_argument("var", "contain non-negative variance estimates")
  }
  df <- rep_len(df, length(var))
  keep <- is.finite(var) & is.finite(df) & df > 0
  if (sum(keep) < 3L) {
    abort_bad_argument("var, df", "include at least three features with finite variance and finite positive df")
  }
  if (sum(var[keep] > 0) < 2L) {
    abort_bad_argument("var", "include at least two positive usable variances to estimate a prior")
  }
  list(var = var[keep], df = df[keep])
}

#' Average ordinary two-sided t-test power over variance and effect priors
#' @param n prospective sample counts per group
#' @param var_prior inverse-gamma prior scale s0 squared
#' @param df_prior prior degrees of freedom, possibly infinite
#' @param delta_mean mean difference among non-null features
#' @param delta_se standard deviation of those differences
#' @param sig_level per-test significance level
#' @noRd
.fdr_prior_power <- function(n, var_prior, df_prior, delta_mean, delta_se, sig_level) {
  if (is.infinite(df_prior)) {
    precision <- 1 / var_prior
  } else {
    nodes <- if (df_prior < 4) 1024L else 256L
    precision <- stats::qgamma(
      (seq_len(nodes) - 0.5) / nodes,
      shape = df_prior / 2,
      rate = df_prior * var_prior / 2
    )
  }
  effects <- if (.is_zero(delta_se)) {
    delta_mean
  } else {
    stats::qnorm((seq_len(512L) - 0.5) / 512L, mean = delta_mean, sd = delta_se)
  }
  effect_sd <- as.numeric(outer(abs(effects), sqrt(precision)))
  vapply(
    n,
    function(size) {
      residual_df <- 2 * size - 2
      critical <- stats::qt(sig_level / 2, df = residual_df, lower.tail = FALSE)
      ncp <- effect_sd * sqrt(size / 2)
      mean(
        stats::pt(critical, df = residual_df, ncp = ncp, lower.tail = FALSE) +
          stats::pt(-critical, df = residual_df, ncp = ncp)
      )
    },
    numeric(1)
  )
}

#' Plan sample size from a pilot variance prior and a planned effect distribution
#'
#' Fit a variance prior to per-feature pilot variances and residual degrees of
#' freedom, then integrate ordinary two-sided t-test power over that prior and
#' the user's planned distribution of non-null effects.
#'
#' @param var real numeric vector of unweighted, unmoderated within-group
#'   variance estimates on the response scale (usually log2 intensity).
#'   Non-finite values are excluded; finite negative variances are errors.
#' @param df real numeric vector of residual degrees of freedom, either length
#'   one (shared by all features) or length(var). Non-finite and non-positive
#'   values are excluded together with their variance estimate.
#' @param delta_mean one finite mean difference among non-null features, on
#'   the same scale as the pilot response; negative or zero means are allowed.
#' @param delta_se one finite, non-negative standard deviation of the true
#'   differences among non-null features. This is not the standard error of an
#'   estimated mean effect. Zero specifies one fixed planned effect.
#' @param max_n one whole number of at least two: the largest sample count
#'   per group to consider.
#' @inheritParams fdr_sig_level
#' @return A list containing `N` (first integer sample count per group reaching
#'   the target, or `NA_integer_`), `power_achieved` (at `N`, or at `max_n` when
#'   the target is not reached), `power_curve` (a tibble of `n` and `power` for
#'   2 through `max_n`), `sig_level`, `var_prior`, `df_prior`,
#'   `n_features_used`, and the original inputs `var`, `df`, `delta_mean`,
#'   `delta_se`, `fdr`, `power`, `pi0`, `max_n`.
#' @details At least three usable features, including two positive variances,
#'   are needed to estimate the prior with
#'   `limma::squeezeVar(var, df, robust = FALSE)`. Zero variances are retained;
#'   any limma warnings about unreliable or offset variances remain visible.
#'   If the estimated prior parameters are \eqn{s_0^2} and \eqn{d_0}, true
#'   variances follow an inverse-gamma distribution with shape \eqn{d_0/2} and
#'   scale \eqn{d_0 s_0^2/2}. True non-null differences independently follow
#'   `Normal(delta_mean, delta_se)` (mean and standard deviation).
#'
#'   Deterministic midpoint quadrature uses gamma-precision and normal-effect
#'   quantiles: 256 variance nodes (1024 when prior df is below four) and 512
#'   effect nodes. Infinite prior df gives the single variance `var_prior`;
#'   zero `delta_se` gives the single effect `delta_mean`.
#'
#'   The power curve uses the fixed threshold from [fdr_sig_level()] at the
#'   requested target power. The chosen sample size achieves the FDR planning
#'   bound under those assumptions; points below the target power need not.
#'   Prospective tests use `2 * n - 2` degrees of freedom, not moderated df.
#'   A moderated test gains prior degrees of freedom, so this ordinary t-test
#'   sample size is conservative relative to that test under the same model.
#'   This function estimates the variance distribution, not effects or `pi0`.
#' @inherit fdr_sig_level references
#' @family stats
#' @export
#' @examples
#' variances <- 0.01 * stats::qf((seq_len(200) - 0.5) / 200, df1 = 8, df2 = 10)
#' planned <- power_fdr_prior_n(variances, df = 8, delta_mean = 0.15, delta_se = 0.02)
#' planned[c("N", "power_achieved", "var_prior", "df_prior", "n_features_used")]
power_fdr_prior_n <- function(var, df, delta_mean, delta_se = 0, fdr = 0.05, power = 0.8, pi0 = 0.95, max_n = 100) {
  sig_level <- fdr_sig_level(fdr = fdr, pi0 = pi0, power = power)
  .check_fdr_effect_inputs(delta_mean, delta_se, max_n)
  usable <- .fdr_prior_inputs(var, df)
  prior <- limma::squeezeVar(usable$var, usable$df, robust = FALSE)
  n <- seq.int(2L, max_n)
  power_curve <- tibble::tibble(
    n = n,
    power = .fdr_prior_power(n, prior$var.prior, prior$df.prior, delta_mean, delta_se, sig_level)
  )
  reached <- match(TRUE, power_curve$power >= power)
  list(
    N = n[reached],
    power_achieved = power_curve$power[if (is.na(reached)) length(n) else reached],
    power_curve = power_curve,
    sig_level = sig_level,
    var_prior = prior$var.prior,
    df_prior = prior$df.prior,
    n_features_used = length(usable$var),
    var = var,
    df = df,
    delta_mean = delta_mean,
    delta_se = delta_se,
    fdr = fdr,
    power = power,
    pi0 = pi0,
    max_n = max_n
  )
}
