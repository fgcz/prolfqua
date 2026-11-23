#' Compute pooled variance
#'
#' following the documentation here:
#' https://online.stat.psu.edu/stat500/lesson/7/7.3/7.3.1/7.3.1.1
#'
#' @param x data.frame
#' @return data.frame
#' @noRd
.compute_pooled <- function(x) {
  xm <- x |> dplyr::filter(.data$nrMeasured > 0)
  mean_all <- sum(xm$meanAbundance * xm$nrMeasured) / sum(xm$nrMeasured)
  x <- x |> dplyr::filter(.data$nrMeasured > 1)
  n <- x$nrMeasured
  pool.n <- sum(n)
  n.groups <- length(x$var)
  pool.var <- sum((n - 1) * x$var) / (pool.n - n.groups)
  pool.mean <- sum(x$meanAbundance * n) / pool.n
  data.frame(
    n.groups = n.groups,
    n = pool.n,
    df = pool.n - n.groups,
    sd = sqrt(pool.var),
    sdT = sqrt(pool.var * 2 / (pool.n / n.groups)),
    var = pool.var,
    mean = if (is.na(pool.mean)) mean_all else pool.mean,
    meanAll = mean_all,
    nrMeasured = sum(xm$nrMeasured)
  )
}

#' pooled variance
#'
#' The values of `.compute_pooled()` for every feature, computed with grouped sums over all features at once.
#' @param res1 per-group statistics from [summarize_stats()]
#' @param config AnalysisConfiguration of the data
#' @return one row per feature, in the order of first appearance in `res1`: the hierarchy keys, the pooled
#'   statistics and the first factor key set to `"pooled"`
#' @noRd
.poolvar <- function(res1, config) {
  grouped <- dplyr::group_by(res1, dplyr::across(dplyr::all_of(config$hierarchy_keys())))
  id <- dplyr::group_indices(grouped)
  first <- id[!duplicated(id)]
  # sums per feature, in row order within each feature (as sum() per feature), features in order of appearance
  sum_by <- function(x) unname(rowsum(x, id, reorder = TRUE)[first, 1])
  measured <- res1$nrMeasured
  observed <- measured > 0L
  replicated <- measured > 1L
  n_measured <- sum_by(measured * observed)
  mean_all <- sum_by(ifelse(observed, res1$meanAbundance * measured, 0)) / n_measured
  pool_n <- sum_by(measured * replicated)
  n_groups <- sum_by(as.integer(replicated))
  pool_var <- sum_by(ifelse(replicated, (measured - 1L) * res1$var, 0)) / (pool_n - n_groups)
  pool_mean <- sum_by(ifelse(replicated, res1$meanAbundance * measured, 0)) / pool_n
  resp <- dplyr::group_keys(grouped)[first, ]
  resp$n.groups <- n_groups
  resp$n <- pool_n
  resp$df <- pool_n - n_groups
  resp$sd <- sqrt(pool_var)
  resp$sdT <- sqrt(pool_var * 2 / (pool_n / n_groups))
  resp$var <- pool_var
  resp$mean <- ifelse(is.na(pool_mean), mean_all, pool_mean)
  resp$meanAll <- mean_all
  resp$nrMeasured <- n_measured
  resp[[config$factor_keys()[1]]] <- "pooled"
  resp
}

#' Counts, sd, var, mean and median of `y` per group, computed for all groups at once
#'
#' Same values as `n()`, `sum(!is.na(y))`, `sd()`, `var()`, `mean()` and `median()` with `na.rm = TRUE` per group:
#' sd and var are NA below two values, the mean is NaN and the median NA for groups without values.
#' @param y numeric response
#' @param group integer group index of each element of `y`, 1 to `n_groups`
#' @param n_groups number of groups
#' @noRd
.group_stats <- function(y, group, n_groups) {
  observed <- !is.na(y)
  y_obs <- y[observed]
  g_obs <- group[observed]
  n <- tabulate(group, n_groups)
  n_measured <- tabulate(g_obs, n_groups)
  group_sum <- function(x) {
    s <- numeric(n_groups)
    s[sort(unique(g_obs))] <- rowsum(x, g_obs, reorder = TRUE)[, 1]
    s
  }
  mean_g <- group_sum(y_obs) / n_measured
  var_g <- group_sum((y_obs - mean_g[g_obs])^2) / (n_measured - 1)
  var_g[n_measured < 2] <- NA_real_
  # medians: values sorted within groups, groups in index order
  sorted <- y_obs[order(g_obs, y_obs)]
  start <- cumsum(c(0L, n_measured))[seq_len(n_groups)]
  median_g <- rep(NA_real_, n_groups)
  has <- which(n_measured > 0)
  lower <- sorted[start[has] + (n_measured[has] + 1L) %/% 2L]
  upper <- sorted[start[has] + n_measured[has] %/% 2L + 1L]
  median_g[has] <- (lower + upper) / 2
  tibble::tibble(
    nrReplicates = n,
    nrMeasured = n_measured,
    nrNAs = n - n_measured,
    sd = sqrt(var_g),
    var = var_g,
    meanAbundance = mean_g,
    medianAbundance = median_g
  )
}

#' Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factor_key character vector - factor columns to group by (default: relevant_factor_keys)
#' @return `summarize_stats()`: a tibble with one row per hierarchy element and level of
#'   `factor_key`, with the columns `nrReplicates`, `nrMeasured`, `nrNAs`, `sd`, `var`, `meanAbundance`,
#'   `medianAbundance`, `CV` (untransformed data only) and `interaction`.
#' @noRd
.summarize_stats <- function(lfqdata, factor_key = lfqdata$relevant_factor_keys()) {
  grouped <- dplyr::group_by(
    lfqdata$data_long(),
    dplyr::across(all_of(c(lfqdata$hierarchy_keys(), factor_key)))
  )
  hierarchy_factor <- dplyr::bind_cols(
    dplyr::group_keys(grouped),
    .group_stats(as.double(grouped[[lfqdata$response()]]), dplyr::group_indices(grouped), dplyr::n_groups(grouped))
  ) |>
    dplyr::mutate(dplyr::across(all_of(factor_key), as.character))
  if (!lfqdata$is_transformed()) {
    hierarchy_factor <- hierarchy_factor |> dplyr::mutate(CV = sd / meanAbundance * 100)
  }
  if (length(factor_key) == 0) {
    hierarchy_factor <- dplyr::mutate(hierarchy_factor, !!lfqdata$factor_keys()[1] := "All")
    hierarchy_factor$interaction <- "All"
  } else {
    hierarchy_factor <- .make_interaction_column(hierarchy_factor, columns = factor_key, sep = ":")
  }
  hierarchy_factor
}


#' Compute variance, sd and related statistics for all factor levels
#'
#' Summary statistics per protein for every level of the relevant factors and of their interaction, in one table;
#' used to judge variability and missingness per group.
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return a tibble with one row per hierarchy element and factor level: the statistics of
#'   `LFQData$get_Stats()$stats()` (`nrReplicates`, `nrMeasured`, `nrNAs`, `sd`, `var`, `meanAbundance`,
#'   `medianAbundance`, `CV` for untransformed data) for the interaction of all relevant factors and, with more than
#'   one factor, for each factor alone; `interaction` names the level.
#' @export
#' @examples
#' res2 <- prolfqua::sim_lfq_data_2factor_config()
#' res2$config$factor_depth <- 2
#' lfq2 <- LFQData$new(res2$data, res2$config)
#' xx <- summarize_stats_factors(lfq2)
#' stopifnot(nrow(xx) == 80)
#' stopifnot(length(unique(xx$interaction)) == (2 + 2 + 2 * 2))
summarize_stats_factors <- function(lfqdata) {
  rfk <- lfqdata$relevant_factor_keys()
  per_factor <- if (length(rfk) > 1) lapply(rfk, function(factor) .summarize_stats(lfqdata, factor_key = factor))
  dplyr::bind_rows(c(list(.summarize_stats(lfqdata)), per_factor))
}


#' summarize stats output (compute quantiles)
#' @param stats_res result of running `summarize_stats`
#' @param factor_keys_depth character vector - factor columns at current depth
#' @param stats summarize either sd or CV
#' @param probs for which quantiles 10, 20 etc.
#' @return `summarize_stats_quantiles()`: a list with `long`, a tibble with one row per factor level and quantile
#'   (columns `probs` and `quantiles`), and `wide`, the same quantiles with one column per factor level.
#' @noRd
.summarize_stats_quantiles <- function(
  stats_res,
  factor_keys_depth,
  stats = c("sd", "CV"),
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
) {
  stats <- match.arg(stats)
  if (!is.numeric(probs) || length(probs) == 0L || anyNA(probs) || any(probs < 0 | probs > 1)) {
    abort_bad_argument("probs", "be a vector of numbers in [0, 1]", not = toString(probs))
  }

  sd_quantile_res2 <- stats_res |>
    dplyr::filter(!is.na(.data[[stats]])) |>
    dplyr::reframe(
      probs = probs,
      quantiles = stats::quantile(.data[[stats]], probs, na.rm = TRUE),
      .by = dplyr::all_of(factor_keys_depth)
    ) |>
    dplyr::arrange(dplyr::pick(dplyr::all_of(factor_keys_depth)))

  xx <- sd_quantile_res2 |> tidyr::unite("interaction", dplyr::all_of(factor_keys_depth))
  wide <- xx |> tidyr::pivot_wider(names_from = "interaction", values_from = "quantiles")
  list(long = sd_quantile_res2, wide = wide)
}


#' Two-sample t-test sample sizes for sd and delta (recycled); each sd is first floored at the smallest sd that
#'
#' min.n samples per group can resolve.
#' @param sd response standard deviations
#' @param delta effect sizes to detect
#' @param power target power
#' @param sig.level test significance level
#' @param min.n minimum sample count per group
#' @noRd
.power_t_test_n <- function(sd, delta, power, sig.level, min.n) {
  if (!is.numeric(delta) || length(delta) == 0L || !all(is.finite(delta)) || any(delta <= 0)) {
    abort_bad_argument("delta", "be a vector of positive finite numbers", not = toString(delta))
  }
  .assert_number(power, "be one number strictly between zero and one", function(v) v > 0 && v < 1)
  .assert_number(sig.level, "be one number strictly between zero and one", function(v) v > 0 && v < 1)
  min_sd <- purrr::map_dbl(delta, function(d) {
    power.t.test(delta = d, n = min.n, sd = NULL, power = power, sig.level = sig.level)$sd
  })
  sdtrimmed <- pmax(sd, min_sd)
  N_exact <- purrr::map2_dbl(sdtrimmed, delta, function(s, d) {
    power.t.test(delta = d, sd = s, power = power, sig.level = sig.level)$n
  })
  list(sdtrimmed = sdtrimmed, N_exact = N_exact, N = ceiling(N_exact))
}

#' Estimate sample sizes from standard deviation quantiles
#' @param quantile_sd output of `summarize_stats_quantiles`
#' @param delta effect size you are interested in
#' @param power of test
#' @param sig.level P-Value
#' @param min.n smallest n to determine
#' @return a tibble with the rows of `quantile_sd` repeated for each `delta`, with the added columns `sdtrimmed`
#'   (the sd floored at the smallest sd resolvable with `min.n` samples), `N_exact`, `N` (rounded up) and `delta`
#'
#' @noRd
.lfq_power_t_test_quantiles <-
  function(quantile_sd, delta = c(0.59, 1, 2), power = 0.8, sig.level = 0.05, min.n = 1.5) {
    res <- lapply(delta, function(d) {
      dplyr::bind_cols(quantile_sd, .power_t_test_n(quantile_sd$quantiles, d, power, sig.level, min.n), delta = d)
    })
    bind_rows(res)
  }

#' Compute theoretical sample sizes from factor level standard deviations
#' @param stats_res data.frame `summarize_stats` output
#' @param delta effect size you are interested in
#' @param power of test
#' @param sig.level P-Value
#' @param min.n smallest n to determine
#' @return a tibble with one row per complete row of `stats_res` and value of `delta`, with the added columns
#'   `delta`, `N_exact` (samples per group needed) and `N` (rounded up)
#'
#' @noRd
.lfq_power_t_test_proteins <- function(stats_res, delta = c(0.59, 1, 2), power = 0.8, sig.level = 0.05, min.n = 1.5) {
  sd_delta <- tidyr::crossing(na.omit(stats_res), delta = delta)
  pw <- .power_t_test_n(sd_delta$sd, sd_delta$delta, power, sig.level, min.n)
  dplyr::bind_cols(sd_delta, pw[c("N_exact", "N")])
}

#' Plot the density or ECDF of sd, mean or CV
#' @param pdata data.frame with statistics
#' @param factor_key character - factor column name for colouring
#' @param stat sd, mean or CV
#' @param ggstat either density or ecdf
#' @return a ggplot object with one density or ECDF line of `stat` per level of `factor_key`
#'
#' @noRd
.plot_stat_density <- function(
  pdata,
  factor_key,
  stat = c("CV", "meanAbundance", "sd"),
  ggstat = c("density", "ecdf")
) {
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  p <- ggplot(pdata, aes(x = .data[[stat]], colour = .data[[factor_key]])) +
    theme_prolfqua(grid = TRUE) +
    # statistics are NA for features with too few values; drop them without a warning
    geom_line(stat = ggstat, na.rm = TRUE)
  p
}
#' Plot the density or ECDF of sd, mean or CV for intensities below and above the median
#' @param pdata data.frame with statistics
#' @param factor_key character - factor column name for faceting
#' @param stat sd, mean or CV
#' @param ggstat either density or ecdf
#' @return a ggplot object with density or ECDF lines of `stat` for features above and below the median mean
#'   intensity, faceted by `factor_key`
#'
#' @noRd
.plot_stat_density_median <- function(
  pdata,
  factor_key,
  stat = c("CV", "meanAbundance", "sd"),
  ggstat = c("density", "ecdf")
) {
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  pdata <- pdata |> dplyr::filter(!is.na(.data[[stat]]))
  top50 <- pdata |>
    dplyr::mutate(top = ifelse(meanAbundance > median(meanAbundance, na.rm = TRUE), "top 50", "bottom 50"))
  p <- ggplot(top50, aes(x = .data[[stat]], colour = .data$top)) +
    theme_prolfqua(grid = TRUE) +
    geom_line(stat = ggstat, na.rm = TRUE) +
    facet_wrap(factor_key)
  p
}

#' Violin plot of sd, CV or mean
#'
#' @param pdata data.frame with statistics
#' @param factor_keys_depth character vector - factor columns for grouping
#' @param stat either CV, mean or sd
#' @return a ggplot object with one violin of `stat` per combination of factor levels, with the median marked
#' @noRd
.plot_stat_violin <- function(pdata, factor_keys_depth, stat = c("CV", "meanAbundance", "sd")) {
  stat <- match.arg(stat)
  pdata <- pdata |> tidyr::unite("groups", dplyr::all_of(factor_keys_depth))
  p <- ggplot(pdata, aes(x = .data$groups, y = .data[[stat]])) +
    theme_prolfqua() +
    # statistics are NA for features with too few values; drop them without a warning
    geom_violin(linewidth = 0.3, na.rm = TRUE) +
    ggplot2::stat_summary(fun = median, geom = "point", size = 1, color = "black", na.rm = TRUE)
  p
}
#' Violin plot of sd, CV or mean for intensities below and above the median
#' @param pdata data.frame with statistics
#' @param factor_key character - factor column name for x-axis
#' @param stat either CV, mean or sd
#' @return a ggplot object with violins of `stat` per level of `factor_key` and the median and interquartile range
#'   marked, faceted into features below and above the median mean intensity
#'
#' @noRd
.plot_stat_violin_median <- function(pdata, factor_key, stat = c("CV", "meanAbundance", "sd")) {
  stat <- match.arg(stat)
  pdata <- pdata |> dplyr::filter(!is.na(.data[[stat]]))

  top50 <- pdata |>
    dplyr::mutate(top = ifelse(meanAbundance > median(meanAbundance, na.rm = TRUE), "top 50", "bottom 50"))

  p <- ggplot(top50, aes(x = .data[[factor_key]], y = .data[[stat]])) +
    theme_prolfqua() +
    geom_violin(linewidth = 0.3, na.rm = TRUE) +
    stat_summary(
      fun = stats::median,
      fun.min = function(x) stats::quantile(x, 0.25),
      fun.max = function(x) stats::quantile(x, 0.75),
      geom = "pointrange",
      na.rm = TRUE
    ) +
    facet_wrap("top")
  p
}

#' Plot standard deviation against mean to assess variance stability
#' @param pdata data.frame with statistics
#' @param factor_keys_depth character vector - factor columns for faceting
#' @param size how many points to sample (since scatter plot to slow for all)
#' @return a ggplot scatter plot of `sd` against `meanAbundance` with a loess fit, one facet per combination of
#'   factor levels
#'
#' @noRd
.plot_stdv_vs_mean <- function(pdata, factor_keys_depth, size = 2000) {
  summary <- pdata |>
    group_by(across(all_of(factor_keys_depth))) |>
    dplyr::summarize(n = n(), .groups = "drop")
  size <- min(size, min(summary$n))

  pdata <- pdata |>
    dplyr::group_by(across(all_of(factor_keys_depth))) |>
    dplyr::slice_sample(n = size) |>
    dplyr::ungroup()

  p <- ggplot(pdata, aes(x = meanAbundance, y = sd)) +
    theme_prolfqua(grid = TRUE) +
    # sd is NA for features with fewer than two values, and 0 or NA leave a log scale
    geom_point(na.rm = TRUE) +
    geom_smooth(method = "loess", formula = y ~ x, na.rm = TRUE) +
    facet_wrap(factor_keys_depth) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
}
