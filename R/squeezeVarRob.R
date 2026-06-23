# Code from https://raw.githubusercontent.com/statOmics/MSqRob/master/R/squeezeVarRob.R

trigammaInverse <- function(x) {
  if (!is.numeric(x)) {
    stop("Non-numeric argument to mathematical function")
  }
  if (length(x) == 0) {
    return(numeric(0))
  }
  omit <- is.na(x)
  if (any(omit)) {
    y <- x
    if (any(!omit)) {
      y[!omit] <- Recall(x[!omit])
    }
    return(y)
  }
  omit <- (x < 0)
  if (any(omit)) {
    y <- x
    y[omit] <- NaN
    warning("NaNs produced")
    if (any(!omit)) {
      y[!omit] <- Recall(x[!omit])
    }
    return(y)
  }
  omit <- (x > 1e+07)
  if (any(omit)) {
    y <- x
    y[omit] <- 1 / sqrt(x[omit])
    if (any(!omit)) {
      y[!omit] <- Recall(x[!omit])
    }
    return(y)
  }
  omit <- (x < 1e-06)
  if (any(omit)) {
    y <- x
    y[omit] <- 1 / x[omit]
    if (any(!omit)) {
      y[!omit] <- Recall(x[!omit])
    }
    return(y)
  }
  y <- 0.5 + 1 / x
  iter <- 0
  repeat {
    iter <- iter + 1
    tri <- trigamma(y)
    dif <- tri * (1 - tri / x) / psigamma(y, deriv = 2)
    y <- y + dif
    if (max(-dif / y) < 1e-08) {
      break
    }
    if (iter > 50) {
      warning("Iteration limit exceeded")
      break
    }
  }
  y
}

# Moment estimation of the parameters of a scaled F-distribution.
# The numerator degrees of freedom are given; the denominator is estimated.
# Gordon Smyth and Belinda Phipson, 8 Sept 2002; last revised 25 Jan 2017.
fitFDist_LG <- function(x, df1, covariate = NULL, min_df = 1) {
  n <- length(x)
  if (n <= 1L) {
    return(list(scale = NA, df2 = NA))
  }

  df1_setup <- .fit_f_dist_prepare_df1(df1, n, min_df)
  if (!is.null(df1_setup$early)) {
    return(df1_setup$early)
  }

  covariate_setup <- .fit_f_dist_prepare_covariate(covariate, n, x, df1, min_df)
  if (!is.null(covariate_setup$early)) {
    return(covariate_setup$early)
  }

  filtered <- .fit_f_dist_filter_inputs(x, df1, covariate_setup$covariate, df1_setup$ok)

  if (filtered$nok <= covariate_setup$spline_df) {
    s20 <- NA
    if (!is.null(covariate)) {
      s20 <- rep_len(s20, n)
    }
    return(list(scale = s20, df2 = NA))
  }

  x <- .fit_f_dist_offset_zero_variances(filtered$x)
  moments <- .fit_f_dist_log_moments(
    x,
    filtered$df1,
    filtered$covariate,
    filtered$covariate_not_ok,
    filtered$ok,
    filtered$not_all_ok,
    n,
    covariate_setup$spline_df
  )

  .fit_f_dist_scale_from_moments(
    x,
    filtered$df1,
    moments$emean,
    moments$evar,
    has_covariate = !is.null(covariate)
  )
}

.fit_f_dist_prepare_df1 <- function(df1, n, min_df) {
  ok <- is.finite(df1) & df1 > min_df #1e-15
  if (length(df1) == 1L) {
    if (!ok) {
      return(list(ok = ok, early = list(scale = NA, df2 = NA)))
    }
    ok <- rep_len(TRUE, n)
  } else {
    if (length(df1) != n) stop("x and df1 have different lengths")
  }

  list(ok = ok, early = NULL)
}

.fit_f_dist_prepare_covariate <- function(covariate, n, x, df1, min_df) {
  if (is.null(covariate)) {
    return(list(covariate = NULL, spline_df = 1L, early = NULL))
  }

  if (length(covariate) != n) {
    stop("x and covariate must be of same length")
  }
  if (anyNA(covariate)) {
    stop("NA covariate values not allowed")
  }

  isfin <- is.finite(covariate)
  if (!all(isfin)) {
    if (any(isfin)) {
      r <- range(covariate[isfin])
      covariate[covariate == -Inf] <- r[1] - 1
      covariate[covariate == Inf] <- r[2] + 1
    } else {
      covariate <- sign(covariate)
    }
  }

  spline_df <- min(4L, length(unique(covariate)))
  # If covariate takes only one value, recall with NULL covariate.
  if (spline_df < 2L) {
    out <- fitFDist_LG(x = x, df1 = df1)
    out$scale <- rep_len(out$scale, n)
    return(list(covariate = covariate, spline_df = spline_df, early = out))
  }

  list(covariate = covariate, spline_df = spline_df, early = NULL)
}

.fit_f_dist_filter_inputs <- function(x, df1, covariate, ok) {
  ok <- ok & is.finite(x) & (x > -1e-15)
  not_all_ok <- !all(ok)
  covariate_not_ok <- NULL

  if (not_all_ok) {
    x <- x[ok]
    if (length(df1) > 1L) {
      df1 <- df1[ok]
    }
    if (!is.null(covariate)) {
      covariate_not_ok <- covariate[!ok]
      covariate <- covariate[ok]
    }
  }

  list(
    x = x,
    df1 = df1,
    covariate = covariate,
    covariate_not_ok = covariate_not_ok,
    ok = ok,
    nok = sum(ok),
    not_all_ok = not_all_ok
  )
}

.fit_f_dist_offset_zero_variances <- function(x) {
  x <- pmax(x, 0)
  m <- median(x)
  if (m == 0) {
    warning("More than half of residual variances are exactly zero: eBayes unreliable")
    m <- 1
  } else {
    if (any(x == 0)) warning("Zero sample variances detected, have been offset", call. = FALSE)
  }
  pmax(x, 1e-5 * m)
}

.fit_f_dist_log_moments <- function(x, df1, covariate, covariate_not_ok, ok, not_all_ok, n, spline_df) {
  z <- log(x)
  e <- z - digamma(df1 / 2) + log(df1 / 2)

  if (is.null(covariate)) {
    emean <- mean(e)
    evar <- sum((e - emean)^2) / (length(x) - 1L) # equals evar <- var(e) if emean=mean(e)
  } else {
    if (!requireNamespace("splines", quietly = TRUE)) {
      stop("splines package required but is not available")
    }
    design <- try(splines::ns(covariate, df = spline_df, intercept = TRUE), silent = TRUE)
    if (inherits(design, "try-error")) {
      stop("Problem with covariate")
    }
    fit <- lm.fit(design, e)
    if (not_all_ok) {
      design2 <- predict(design, newx = covariate_not_ok)
      emean <- rep_len(0, n)
      emean[ok] <- fit$fitted
      emean[!ok] <- design2 %*% fit$coefficients
    } else {
      emean <- fit$fitted
    }
    evar <- mean(fit$effects[-seq_len(fit$rank)]^2)
  }

  #MSqRob added: avoid NaN in evar
  if (n == 1) {
    evar <- 0
  }

  list(emean = emean, evar = evar)
}

.fit_f_dist_scale_from_moments <- function(x, df1, emean, evar, has_covariate) {
  evar <- evar - mean(trigamma(df1 / 2))
  if (evar > 0) {
    df2 <- 2 * trigammaInverse(evar)
    s20 <- exp(emean + digamma(df2 / 2) - log(df2 / 2))
  } else {
    df2 <- Inf
    if (!has_covariate) {
      #      Use simple pooled variance, which is MLE of the scale in this case.
      #      Versions of limma before Jan 2017 returned the limiting value of the evar>0 estimate, which is larger.
      s20 <- mean(x)
    } else {
      s20 <- exp(emean)
    }
  }

  list(scale = s20, df2 = df2)
}


.validate_fit_f_dist_robustly <- function(x, df1, covariate) {
  n <- length(x)
  if (all(length(df1) != c(1, n))) {
    stop("x and df1 are different lengths")
  }
  if (!is.null(covariate)) {
    if (length(covariate) != n) {
      stop("x and covariate are different lengths")
    }
    if (!all(is.finite(covariate))) {
      stop("covariate contains NA or infinite values")
    }
  }
}

.filter_fit_f_dist_robustly <- function(x, df1, covariate, ok, winsor.tail.p, trace) {
  df2.shrunk <- x
  x_ok <- x[ok]
  df1_ok <- df1
  if (length(df1) > 1) {
    df1_ok <- df1[ok]
  }
  covariate_ok <- covariate
  covariate_not_ok <- NULL
  if (!is.null(covariate)) {
    covariate_not_ok <- covariate[!ok]
    covariate_ok <- covariate[ok]
  }

  fit <- fitFDistRobustly_LG(
    x = x_ok,
    df1 = df1_ok,
    covariate = covariate_ok,
    winsor.tail.p = winsor.tail.p,
    trace = trace
  )
  df2.shrunk[ok] <- fit$df2.shrunk
  df2.shrunk[!ok] <- fit$df2
  if (is.null(covariate)) {
    scale <- fit$scale
  } else {
    scale <- x_ok
    scale[ok] <- fit$scale
    scale[!ok] <- exp(approx(covariate_ok, log(fit$scale), xout = covariate_not_ok, rule = 2)$y)
  }
  list(scale = scale, df2 = fit$df2, df2.shrunk = df2.shrunk)
}

.offset_small_variances <- function(x) {
  m <- median(x)
  if (m <= 0) {
    stop("x values are mostly <= 0")
  }
  i <- (x < m * 1e-12)
  if (any(i)) {
    nzero <- sum(i)
    if (nzero == 1) {
      warning("One very small variance detected, has been offset away from zero", call. = FALSE)
    } else {
      warning(nzero, " very small variances detected, have been offset away from zero", call. = FALSE)
    }
    x[i] <- m * 1e-12
  }
  x
}

.transform_to_constant_df1 <- function(x, df1, covariate, non_robust) {
  if (length(df1) <= 1) {
    return(list(x = x, df1 = df1))
  }

  df1max <- max(df1)
  i <- (df1 < (df1max - 1e-14))
  if (!any(i)) {
    return(list(x = x, df1 = df1[1]))
  }

  if (is.null(covariate)) {
    s <- non_robust$scale
  } else {
    s <- non_robust$scale[i]
  }
  f <- x[i] / s
  df2 <- non_robust$df2
  pupper <- pf(f, df1 = df1[i], df2 = df2, lower.tail = FALSE, log.p = TRUE)
  plower <- pf(f, df1 = df1[i], df2 = df2, lower.tail = TRUE, log.p = TRUE)
  up <- pupper < plower
  if (any(up)) {
    f[up] <- qf(pupper[up], df1 = df1max, df2 = df2, lower.tail = FALSE, log.p = TRUE)
  }
  if (any(!up)) {
    f[!up] <- qf(plower[!up], df1 = df1max, df2 = df2, lower.tail = TRUE, log.p = TRUE)
  }
  x[i] <- f * s
  list(x = x, df1 = df1max)
}

.winsorized_residual_moments <- function(z, covariate, winsor.tail.p, prob, trace) {
  if (is.null(covariate)) {
    ztrend <- mean(z, trim = winsor.tail.p[2])
    zresid <- z - ztrend
  } else {
    lo <- limma::loessFit(z, covariate, span = 0.4)
    ztrend <- lo$fitted
    zresid <- lo$residual
  }

  zrq <- quantile(zresid, prob = prob)
  zwins <- pmin(pmax(zresid, zrq[1]), zrq[2])
  zwmean <- median(zwins)
  zwvar <- mad(zwins)
  if (trace) {
    message("Variance of Winsorized Fisher-z ", zwvar)
  }
  list(ztrend = ztrend, zwmean = zwmean, zwvar = zwvar)
}

.winsorized_moments_function <- function() {
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("statmod package required but is not installed")
  }
  g <- statmod::gauss.quad.prob(128, dist = "uniform")
  linkfun <- function(x) x / (1 + x)
  linkinv <- function(x) x / (1 - x)

  list(
    linkfun = linkfun,
    linkinv = linkinv,
    moments = function(df1, df2, winsor.tail.p) {
      fq <- qf(p = c(winsor.tail.p[1], 1 - winsor.tail.p[2]), df1 = df1, df2 = df2)
      zq <- log(fq)
      q <- linkfun(fq)
      nodes <- q[1] + (q[2] - q[1]) * g$nodes
      fnodes <- linkinv(nodes)
      znodes <- log(fnodes)
      f <- df(fnodes, df1 = df1, df2 = df2) / (1 - nodes)^2
      q21 <- q[2] - q[1]
      m <- q21 * sum(g$weights * f * znodes) + sum(zq * winsor.tail.p)
      v <- q21 * sum(g$weights * f * (znodes - m)^2) + sum((zq - m)^2 * winsor.tail.p)
      list(mean = m, var = v)
    }
  )
}

.fit_infinite_df2 <- function(z, ztrend, zwmean, mom, df1, n) {
  df2 <- Inf
  ztrendcorrected <- ztrend + zwmean - mom$mean
  s20 <- exp(ztrendcorrected)
  fstat <- exp(z - ztrendcorrected)
  tail_p <- pchisq(fstat * df1, df = df1, lower.tail = FALSE)
  r <- rank(fstat)
  empirical_tail_prob <- (n - r + 0.5) / n
  prob_not_outlier <- pmin(tail_p / empirical_tail_prob, 1)
  df_pooled <- n * df1
  df2.shrunk <- rep.int(df2, n)
  outlier <- prob_not_outlier < 1
  if (any(outlier)) {
    df2.shrunk[outlier] <- prob_not_outlier[outlier] * df_pooled
    o <- order(tail_p)
    df2.shrunk[o] <- cummax(df2.shrunk[o])
  }
  list(scale = s20, df2 = df2, df2.shrunk = df2.shrunk)
}

.estimate_robust_df2 <- function(non_robust, df1, winsor.tail.p, zwvar, funval_inf, moments, trace, n) {
  if (non_robust$df2 == Inf) {
    non_robust$df2.shrunk <- rep.int(non_robust$df2, n)
    return(non_robust)
  }

  fun <- function(x) {
    df2 <- moments$linkinv(x)
    mom <- moments$moments(df1 = df1, df2 = df2, winsor.tail.p = winsor.tail.p)
    if (trace) {
      message("df2=", df2, ", Working Var=", mom$var)
    }
    log(zwvar / mom$var)
  }

  rbx <- moments$linkfun(non_robust$df2)
  funval_low <- fun(rbx)
  if (funval_low >= 0) {
    return(non_robust$df2)
  }

  u <- uniroot(fun, interval = c(rbx, 1), tol = 1e-8, f.lower = funval_low, f.upper = funval_inf)
  moments$linkinv(u$root)
}

.shrink_outlier_df2 <- function(z, ztrendcorrected, df1, df2, n) {
  zresid <- z - ztrendcorrected
  fstat <- exp(zresid)
  log_tail_p <- pf(fstat, df1 = df1, df2 = df2, lower.tail = FALSE, log.p = TRUE)
  tail_p <- exp(log_tail_p)
  r <- rank(fstat)
  log_empirical_tail_prob <- log(n - r + 0.5) - log(n)
  log_prob_not_outlier <- pmin(log_tail_p - log_empirical_tail_prob, 0)
  prob_not_outlier <- exp(log_prob_not_outlier)
  prob_outlier <- -expm1(log_prob_not_outlier)

  if (any(prob_not_outlier < 1)) {
    min_log_tail_p <- min(log_tail_p)
    if (min_log_tail_p == -Inf) {
      df2.outlier <- 0
      df2.shrunk <- prob_not_outlier * df2
    } else {
      df2.outlier <- log(0.5) / min_log_tail_p * df2
      new_log_tail_p <- pf(max(fstat), df1 = df1, df2 = df2.outlier, lower.tail = FALSE, log.p = TRUE)
      df2.outlier <- log(0.5) / new_log_tail_p * df2.outlier
      df2.shrunk <- prob_not_outlier * df2 + prob_outlier * df2.outlier
    }

    o <- order(log_tail_p)
    df2.ordered <- df2.shrunk[o]
    m <- cumsum(df2.ordered)
    m <- m / seq_len(n)
    imin <- which.min(m)
    df2.ordered[seq_len(imin)] <- m[imin]
    df2.shrunk[o] <- cummax(df2.ordered)
  } else {
    df2.outlier <- df2
    df2.shrunk <- rep.int(df2, n)
  }

  list(
    tail.p.value = tail_p,
    prob.outlier = prob_outlier,
    df2.outlier = df2.outlier,
    df2.shrunk = df2.shrunk
  )
}

.fit_f_dist_robustly_complete <- function(x, df1, covariate, winsor.tail.p, trace, min_df, n) {
  x <- .offset_small_variances(x)
  non_robust <- fitFDist_LG(x = x, df1 = df1, covariate = covariate, min_df = min_df)

  prob <- winsor.tail.p <- rep(winsor.tail.p, length = 2)
  prob[2] <- 1 - winsor.tail.p[2]
  if (all(winsor.tail.p < 1 / n)) {
    non_robust$df2.shrunk <- rep.int(non_robust$df2, n)
    return(non_robust)
  }

  transformed <- .transform_to_constant_df1(x, df1, covariate, non_robust)
  x <- transformed$x
  df1 <- transformed$df1

  z <- log(x)
  residual_moments <- .winsorized_residual_moments(z, covariate, winsor.tail.p, prob, trace)
  moments <- .winsorized_moments_function()
  mom <- moments$moments(df1 = df1, df2 = Inf, winsor.tail.p = winsor.tail.p)
  funval_inf <- log(residual_moments$zwvar / mom$var)
  if (funval_inf <= 0) {
    return(.fit_infinite_df2(z, residual_moments$ztrend, residual_moments$zwmean, mom, df1, n))
  }

  df2 <- .estimate_robust_df2(
    non_robust,
    df1,
    winsor.tail.p,
    residual_moments$zwvar,
    funval_inf,
    moments,
    trace,
    n
  )
  if (is.list(df2)) {
    return(df2)
  }

  mom <- moments$moments(df1 = df1, df2 = df2, winsor.tail.p = winsor.tail.p)
  ztrendcorrected <- residual_moments$ztrend + residual_moments$zwmean - mom$mean
  outlier_fit <- .shrink_outlier_df2(z, ztrendcorrected, df1, df2, n)
  c(list(scale = exp(ztrendcorrected), df2 = df2), outlier_fit)
}

# Just adapted to use fitFDist_LG internally instead of fitFDist
fitFDistRobustly_LG <- function(x, df1, covariate = NULL, winsor.tail.p = c(0.05, 0.1), trace = FALSE, min_df = 1) {
  n <- length(x)

  if (n < 2) {
    return(list(scale = NA, df2 = NA, df2.shrunk = NA))
  }
  if (n == 2) {
    return(fitFDist_LG(x = x, df1 = df1, covariate = covariate, min_df = min_df))
  }

  .validate_fit_f_dist_robustly(x, df1, covariate)

  ok <- !is.na(x) & is.finite(df1) & (df1 > min_df)
  if (!all(ok)) {
    return(.filter_fit_f_dist_robustly(x, df1, covariate, ok, winsor.tail.p, trace))
  }

  .fit_f_dist_robustly_complete(x, df1, covariate, winsor.tail.p, trace, min_df, n)
}

.squeezeVarRob <- function(var, df, var.prior, df.prior) {
  #  Squeeze posterior variances given hyperparameters
  #  NAs not allowed in df.prior
  #  Gordon Smyth
  #  Created 5 May 2016

  n <- length(var)
  isfin <- is.finite(df.prior)
  if (all(isfin)) {
    return((df * var + df.prior * var.prior) / (df + df.prior))
  }

  #  From here, at least some df.prior are infinite

  #  For infinite df.prior, return var.prior
  if (length(var.prior) == n) {
    var.post <- var.prior
  } else {
    var.post <- rep_len(var.prior, length.out = n)
  }

  #  Maybe some df.prior are finite
  if (any(isfin)) {
    i <- which(isfin)
    if (length(df) > 1) {
      df <- df[i]
    }
    df.prior <- df.prior[i]
    var.post[i] <- (df * var[i] + df.prior * var.post[i]) / (df + df.prior)
  }

  return(var.post)
}


#' Robustly Squeeze Sample Variances
#'
#' @description This function squeezes a set of sample variances together by
#'   computing empirical Bayes posterior means in a way that is robust against
#'   the presence of very small non-integer degrees of freedom values.
#' @param var A numeric vector of independent sample variances.
#' @param df A numeric vector of degrees of freedom for the sample variances.
#' @param covariate If \code{non-NULL}, \code{var.prior} will depend on this
#'   numeric covariate. Otherwise, \code{var.prior} is constant.
#' @param robust A logical indicating wheter the estimation of
#'   \code{df.prior} and \code{var.prior} should be robustified against outlier
#'   sample variances. Defaults to \code{FALSE}.
#' @param winsor.tail.p A numeric vector of length 1 or 2, giving left and
#'   right tail proportions of \code{x} to Winsorize. Used only when
#'   \code{robust=TRUE}.
#' @param min_df A numeric value indicating the minimal degrees of freedom that
#'   will be taken into account for calculating the prior degrees of freedom
#'   and prior variance.
#' @return A list with components:
#' \code{var.post} A numeric vector of posterior variances.
#' \code{var.prior} The location of prior distribution. A vector if
#'   \code{covariate} is non-\code{NULL}, otherwise a scalar.
#' \code{df.prior} The degrees of freedom of prior distribution. A vector if \code{robust=TRUE}, otherwise a scalar.
#' @export
#' @examples
#' var <- rexp(1000)
#' df <- sample( 3:10, 1000, replace = TRUE)
#' tmp <- squeezeVarRob(var, df)
#' tmp <- squeezeVarRob(var, df, robust = TRUE)
#'
squeezeVarRob <- function(var, df, covariate = NULL, robust = FALSE, winsor.tail.p = c(0.05, 0.1), min_df = 1) {
  # Empirical Bayes posterior variances
  #  Adapted from Gordon Smyth
  # Based on version created 2 March 2004.  Last modified 5 May 2016.
  n <- length(var)

  #  Degenerate special case
  if (n == 0) {
    stop("var is empty")
  }
  # Historical special case removed by MSqRob; keep identical results when
  # there is only one non-missing observation.

  #If there is only one variance that is not NA: we want the same results if we have only one observation!
  if (sum(!is.na(var)) == 1) {
    return(list(var.post = var, var.prior = var, df.prior = df))
  }

  # Keep NA at NA and Inf at Inf (no guard against missing/infinite values when df==0)

  # If only one df is given, repeat it
  if (length(df) == 1) {
    df <- rep.int(df, n)
  } else {
    if (length(df) != n) stop("lengths differ")
  }

  #  Estimate prior var and df
  if (robust) {
    fit <- fitFDistRobustly_LG(var, df1 = df, covariate = covariate, winsor.tail.p = winsor.tail.p, min_df = min_df)
    df.prior <- fit$df2.shrunk
  } else {
    fit <- fitFDist_LG(var, df1 = df, covariate = covariate, min_df = min_df)
    df.prior <- fit$df2
  }
  #  if(anyNA(df.prior)) stop("Could not estimate prior df") -> we want NA's to stay where they are!

  #  Posterior variances
  var.post <- .squeezeVarRob(var = var, df = df, var.prior = fit$scale, df.prior = df.prior)

  list(df.prior = df.prior, var.prior = fit$scale, var.post = var.post)
}
