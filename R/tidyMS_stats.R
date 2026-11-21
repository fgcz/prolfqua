#' Compute pooled variance
#'
#' following the documentation here:
#' https://online.stat.psu.edu/stat500/lesson/7/7.3/7.3.1/7.3.1.1
#'
#' @param x data.frame
#' @return data.frame
#' @export
#' @rdname pooled_var
#' @keywords internal
#' @family stats
#'
#' @examples
#' x <- data.frame(nrMeasured =c(1,2,2), var = c(3,4,4), meanAbundance = c(3,3,3))
#' x <- data.frame(nrMeasured = c(1,2,1,1), var = c(NA, 0.0370, NA, NA),
#'   meanAbundance = c(-1.94,-1.46,-1.87,-1.45))
#' compute_pooled(x)
#' y <- data.frame(dilution.=c("a","b","c"),
#'      nrReplicates = c(4,4,4), nrMeasured = c(0,0,1), sd =c(NA,NA,NA),
#'      var = c(NA,NA,NA),meanAbundance = c(NaN,NaN,NaN))
#' compute_pooled(y)
#' yb <- y |> dplyr::filter(nrMeasured > 1)
compute_pooled <- function(x) {
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
#' @export
#' @rdname pooled_var
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb$data, bb$config)
#' res1 <- summarize_stats(lfq)
#' pv <- poolvar(res1, bb$config)
#' stopifnot(nrow(pv) == nrow(res1) / 3)
#'
poolvar <- function(res1, config) {
  resp <- res1 |> nest(data = -all_of(config$hierarchy_keys()))
  pooled <- purrr::list_rbind(purrr::map(resp$data, compute_pooled))
  resp$data <- NULL
  resp <- bind_cols(resp, pooled)
  resp <- resp |> mutate(!!config$factor_keys()[1] := "pooled")
  resp
}

#' Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @param factor_key character vector - factor columns to group by (default: relevant_factor_keys)
#' @return `summarize_stats()`: a tibble with one row per hierarchy element, isotope label and level of
#'   `factor_key`, with the columns `nrReplicates`, `nrMeasured`, `nrNAs`, `sd`, `var`, `meanAbundance`,
#'   `medianAbundance`, `CV` (untransformed data only) and `interaction`.
#' @export
#' @rdname summarize_stats
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_protein_config()
#' lfq <- LFQData$new(bb$data, bb$config)
#' res1 <- summarize_stats(lfq)
#'
#' res2 <- prolfqua::sim_lfq_data_2factor_config()
#' res2$config$factor_depth <- 2
#' lfq2 <- LFQData$new(res2$data, res2$config)
#' stats <- summarize_stats(lfq2)
#' stopifnot(nrow(stats) == 40)
#'
#' stats <- summarize_stats(lfq2, factor_key = lfq2$factor_keys()[1])
#' stopifnot(nrow(stats) == 20)
#' stats <- summarize_stats(lfq2, factor_key = lfq2$factor_keys()[2])
#' stopifnot(nrow(stats) == 20)
#' stats <- summarize_stats(lfq2, factor_key = NULL)
#' stopifnot(nrow(stats) == 10)
summarize_stats <- function(lfqdata, factor_key = lfqdata$relevant_factor_keys()) {
  intsym <- sym(lfqdata$response())
  hierarchy_factor <- lfqdata$data_long() |>
    dplyr::group_by(!!!syms(c(lfqdata$hierarchy_keys(), lfqdata$isotope_label(), factor_key))) |>
    dplyr::summarize(
      nrReplicates = dplyr::n(),
      nrMeasured = sum(!is.na(!!intsym)),
      nrNAs = sum(is.na(!!intsym)),
      sd = stats::sd(!!intsym, na.rm = TRUE),
      var = stats::var(!!intsym, na.rm = TRUE),
      meanAbundance = mean(!!intsym, na.rm = TRUE),
      medianAbundance = median(!!intsym, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(all_of(factor_key), as.character))
  if (!lfqdata$is_transformed()) {
    hierarchy_factor <- hierarchy_factor |> dplyr::mutate(CV = sd / meanAbundance * 100)
  }
  if (length(factor_key) == 0) {
    hierarchy_factor <- dplyr::mutate(hierarchy_factor, !!lfqdata$factor_keys()[1] := "All")
    hierarchy_factor$interaction <- "All"
  } else {
    hierarchy_factor <- prolfqua::make_interaction_column(hierarchy_factor, columns = factor_key, sep = ":")
  }
  hierarchy_factor
}


#' Compute variance, sd and related statistics for all factor levels
#'
#' @param lfqdata [LFQData] object with the data and its configuration
#' @return The computed result.
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
  per_factor <- if (length(rfk) > 1) lapply(rfk, function(factor) summarize_stats(lfqdata, factor_key = factor))
  dplyr::bind_rows(c(list(summarize_stats(lfqdata)), per_factor))
}


#' summarize stats output (compute quantiles)
#' @param stats_res result of running `summarize_stats`
#' @param factor_keys_depth character vector - factor columns at current depth
#' @param stats summarize either sd or CV
#' @param probs for which quantiles 10, 20 etc.
#' @return `summarize_stats_quantiles()`: a list with `long`, a tibble with one row per factor level and quantile
#'   (columns `probs` and `quantiles`), and `wide`, the same quantiles with one column per factor level.
#' @rdname summarize_stats
#' @export
#' @keywords internal
#' @family stats
#' @examples
#' library(ggplot2)
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb1$data, bb1$config)
#' stats_res <- summarize_stats(lfq)
#' sq <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys())
#' sq <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys(), stats = "CV")
#' sq <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys(), stats = "sd")
#' xx <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys(), probs = seq(0, 1, by = 0.1))
#' ggplot2::ggplot(xx$long, aes(x = probs, y = quantiles, color = group_)) + geom_line() + geom_point()
#'
summarize_stats_quantiles <- function(
  stats_res,
  factor_keys_depth,
  stats = c("sd", "CV"),
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
) {
  stats <- match.arg(stats)

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
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#'
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb1$data, bb1$config)
#' stats_res <- summarize_stats(lfq)
#' xx <- summarize_stats_quantiles(stats_res, lfq$relevant_factor_keys(), probs = c(0.5, 0.8))
#' bbb <- lfq_power_t_test_quantiles(xx$long)
#' bbb <- dplyr::bind_rows(bbb)
#' summary <- bbb |>
#'  dplyr::select( -N_exact, -quantiles, -sdtrimmed ) |>
#'  tidyr::pivot_wider(names_from = delta, values_from = N)
#'
lfq_power_t_test_quantiles <-
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
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' ldata <- LFQData$new(bb1$data, bb1$config)
#' stats_res <- summarize_stats(ldata)
#' bb <- lfq_power_t_test_proteins(stats_res)
#'
lfq_power_t_test_proteins <- function(stats_res, delta = c(0.59, 1, 2), power = 0.8, sig.level = 0.05, min.n = 1.5) {
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
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb1$data, bb1$config)
#' res <- lfq$get_Stats()$stats()
#' plot_stat_density(res, lfq$factor_keys()[1], stat = "meanAbundance")
#' plot_stat_density(res, lfq$factor_keys()[1], stat = "sd")
#' plot_stat_density(res, lfq$factor_keys()[1], stat = "CV")
plot_stat_density <- function(pdata, factor_key, stat = c("CV", "meanAbundance", "sd"), ggstat = c("density", "ecdf")) {
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  p <- ggplot(pdata, aes(x = .data[[stat]], colour = .data[[factor_key]])) +
    theme_prolfqua(grid = TRUE) +
    geom_line(stat = ggstat)
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
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb1$data, bb1$config)
#' res <- lfq$get_Stats()$stats()
#' plot_stat_density_median(res, lfq$factor_keys()[1], "CV")
#' plot_stat_density_median(res, lfq$factor_keys()[1], "meanAbundance")
#' plot_stat_density_median(res, lfq$factor_keys()[1], "sd")
plot_stat_density_median <- function(
  pdata,
  factor_key,
  stat = c("CV", "meanAbundance", "sd"),
  ggstat = c("density", "ecdf")
) {
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  pdata <- pdata |> dplyr::filter(!is.na(!!sym(stat)))
  top50 <- pdata |>
    dplyr::mutate(top = ifelse(meanAbundance > median(meanAbundance, na.rm = TRUE), "top 50", "bottom 50"))
  p <- ggplot(top50, aes(x = .data[[stat]], colour = .data$top)) +
    theme_prolfqua(grid = TRUE) +
    geom_line(stat = ggstat) +
    facet_wrap(factor_key)
  p
}

#' Violin plot of sd, CV or mean
#'
#' @param pdata data.frame with statistics
#' @param factor_keys_depth character vector - factor columns for grouping
#' @param stat either CV, mean or sd
#' @return a ggplot object with one violin of `stat` per combination of factor levels, with the median marked
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb1$data, bb1$config)
#' res <- lfq$get_Stats()$stats()
#' plot_stat_violin(res, lfq$relevant_factor_keys(), stat = "meanAbundance")
#' plot_stat_violin(res, lfq$relevant_factor_keys(), stat = "sd")
#' plot_stat_violin(res, lfq$relevant_factor_keys(), stat = "CV")
#'
plot_stat_violin <- function(pdata, factor_keys_depth, stat = c("CV", "meanAbundance", "sd")) {
  stat <- match.arg(stat)
  pdata <- pdata |> tidyr::unite("groups", dplyr::all_of(factor_keys_depth))
  p <- ggplot(pdata, aes(x = .data$groups, y = .data[[stat]])) +
    theme_prolfqua() +
    geom_violin(linewidth = 0.3) +
    ggplot2::stat_summary(fun = median, geom = "point", size = 1, color = "black")
  p
}
#' Violin plot of sd, CV or mean for intensities below and above the median
#' @param pdata data.frame with statistics
#' @param factor_key character - factor column name for x-axis
#' @param stat either CV, mean or sd
#' @return a ggplot object with violins of `stat` per level of `factor_key` and the median and interquartile range
#'   marked, faceted into features below and above the median mean intensity
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb1$data, bb1$config)
#' res <- lfq$get_Stats()$stats()
#' plot_stat_violin_median(res, lfq$factor_keys()[1], stat = "meanAbundance")
plot_stat_violin_median <- function(pdata, factor_key, stat = c("CV", "meanAbundance", "sd")) {
  stat <- match.arg(stat)
  pdata <- pdata |> dplyr::filter(!is.na(!!sym(stat)))

  top50 <- pdata |>
    dplyr::mutate(top = ifelse(meanAbundance > median(meanAbundance, na.rm = TRUE), "top 50", "bottom 50"))

  p <- ggplot(top50, aes(x = .data[[factor_key]], y = .data[[stat]])) +
    theme_prolfqua() +
    geom_violin(linewidth = 0.3) +
    stat_summary(
      fun = stats::median,
      fun.min = function(x) stats::quantile(x, 0.25),
      fun.max = function(x) stats::quantile(x, 0.75),
      geom = "pointrange"
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
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb1$data, bb1$config)
#' res <- lfq$get_Stats()$stats()
#' plot_stdv_vs_mean(res, lfq$relevant_factor_keys())
#'
plot_stdv_vs_mean <- function(pdata, factor_keys_depth, size = 2000) {
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
    geom_point() +
    geom_smooth(method = "loess", formula = y ~ x) +
    facet_wrap(factor_keys_depth) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
}
