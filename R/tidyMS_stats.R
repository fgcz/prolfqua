#' compute pooled variance
#' @rdname pooled_var
#' @param x data.frame
#' @return data.frame
#' @examples
#'
#' x <- data.frame(nrMeasured =c(1,2,2), var = c(3,4,4), meanAbundance = c(3,3,3))
#' x <- data.frame(nrMeasured = c(1,2,1,1), var = c(NA, 0.0370, NA, NA),
#'   meanAbundance = c(-1.94,-1.46,-1.87,-1.45))
#' prolfqua:::pooled_V2(na.omit(x))
#' prolfqua:::pooled_V1(na.omit(x))
#' x <- x[1,, drop=FALSE]
#' x
#' na.omit(x)
#' prolfqua:::pooled_V2(na.omit(x))
pooled_V2 <- function(x) {
  n <- x$nrMeasured
  sample.var <- x$var
  sample.mean <- x$meanAbundance
  pool.n <- sum(n)

  pool.mean <- sum(n * sample.mean) / pool.n
  deviation <- sample.mean - pool.mean

  SS <- (n - 1) * sample.var
  pool.SS <- sum(SS) + sum(n * deviation^2)
  pool.var <- pool.SS / (pool.n - 1)
  n.groups <- length(sample.var)
  sd_total <- sqrt(pool.var * 2 / (pool.n / n.groups))

  res <- data.frame(
    n.groups = n.groups,
    n = pool.n,
    df = pool.n - n.groups,
    sd = sqrt(pool.var),
    var = pool.var,
    sdT = sd_total,
    mean = pool.mean
  )
  return(res)
}

#' compute pooled variance V1
#' @rdname pooled_var
#' @param x data.frame
pooled_V1 <- function(x) {
  n <- x$nrMeasured
  sample.var <- x$var
  sample.mean <- x$meanAbundance
  pool.n <- sum(n)

  n.groups <- length(sample.var)
  SS <- (n - 1) * sample.var
  pool.var <- sum(SS) / (pool.n - n.groups)

  pool.mean <- sum(sample.mean * n) / pool.n

  sd_total <- sqrt(pool.var * 2 / (pool.n / n.groups))

  res <- data.frame(
    n.groups = n.groups,
    n = pool.n,
    df = pool.n - n.groups,
    sd = sqrt(pool.var),
    sdT = sd_total,
    var = pool.var,
    mean = pool.mean
  )
  return(res)
}

#' compute pooled variance
#'
#' following the documentation here:
#' https://online.stat.psu.edu/stat500/lesson/7/7.3/7.3.1/7.3.1.1
#'
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
#' compute_pooled(x, method = "V2")
#' y <- data.frame(dilution.=c("a","b","c"),
#'      nrReplicates = c(4,4,4), nrMeasured = c(0,0,1), sd =c(NA,NA,NA),
#'      var = c(NA,NA,NA),meanAbundance = c(NaN,NaN,NaN))
#' compute_pooled(y)
#' yb <- y |> dplyr::filter(nrMeasured > 1)
compute_pooled <- function(x, method = c("V1", "V2")) {
  method <- match.arg(method)
  xm <- x |> dplyr::filter(.data$nrMeasured > 0)
  mean_all <- sum(xm$meanAbundance * xm$nrMeasured) / sum(xm$nrMeasured)
  nr_measured <- sum(xm$nrMeasured)

  func <- pooled_V1
  if (method == "V2") {
    func <- pooled_V2
  }
  x <- x |> dplyr::filter(.data$nrMeasured > 1)

  res <- func(x)
  if (is.na(res$mean)) {
    res$mean <- mean_all
  }
  res$meanAll <- mean_all
  res$nrMeasured <- nr_measured
  return(res)
}

#' pooled variance
#' @export
#' @rdname pooled_var
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' config <- bb$config
#' data <- bb$data
#'
#' res1 <- summarize_stats(data, config)
#' pv <- poolvar(res1, config)
#' stopifnot(nrow(pv) == nrow(res1)/3)
#'
poolvar <- function(res1, config, method = c("V1", "V2")) {
  method <- match.arg(method)
  resp <- res1 |> nest(data = -all_of(config$hierarchy_keys()))
  pooled <- purrr::map_df(resp$data, compute_pooled, method = method)
  resp$data <- NULL
  resp <- bind_cols(resp, pooled)
  resp <- resp |> mutate(!!config$factor_keys()[1] := "pooled")
  return(resp)
}

#' Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @rdname summarize_stats
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' bb <- prolfqua::sim_lfq_data_protein_config()
#' config <- bb$config
#' data <- bb$data
#'
#' res1 <- summarize_stats(data, config)
#'
#' res2 <- prolfqua::sim_lfq_data_2factor_config()
#' res2$config$factor_depth <- 2
#' stats <- summarize_stats(res2$data, res2$config)
#' stopifnot(nrow(stats) == 40)
#'
#' stats <- summarize_stats(res2$data, res2$config, factor_key = res2$config$factor_keys()[1])
#' stopifnot(nrow(stats) == 20)
#' stats <- summarize_stats(res2$data, res2$config, factor_key = res2$config$factor_keys()[2])
#' stopifnot(nrow(stats) == 20)
#' stats <- summarize_stats(res2$data, res2$config, factor_key = NULL)
#' stopifnot(nrow(stats) == 10)
#' # TODO (WEW) add test when there is one level per group.
summarize_stats <- function(pdata, config, factor_key = config$factor_keys_depth(), .completed = FALSE) {
  if (!.completed) {
    pdata <- complete_cases(pdata, config)
  }
  intsym <- sym(config$get_response())
  hierarchy_factor <- pdata |>
    dplyr::group_by(!!!syms(c(config$hierarchy_keys(), config$isotope_label, factor_key))) |>
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
    dplyr::ungroup()

  hierarchy_factor <- hierarchy_factor |>
    dplyr::mutate(dplyr::across(all_of(factor_key), as.character))
  if (config$is_response_transformed == FALSE) {
    hierarchy_factor <- hierarchy_factor |> dplyr::mutate(CV = sd / meanAbundance * 100)
  }
  if (is.null(factor_key) || length(factor_key) == 0) {
    hierarchy_factor <- dplyr::mutate(hierarchy_factor, !!config$factor_keys()[1] := "All")
  }
  hierarchy_factor <- ungroup(hierarchy_factor)
  if (length(factor_key) > 0 && !is.null(factor_key)) {
    hierarchy_factor <- prolfqua::make_interaction_column(hierarchy_factor, columns = factor_key, sep = ":")
  } else {
    hierarchy_factor$interaction <- "All"
  }
  return(hierarchy_factor)
}


#' compute var sd etc for all factor levels
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @examples
#' # example code
#' res2 <- prolfqua::sim_lfq_data_2factor_config()
#' xx <- summarize_stats_factors(res2$data, res2$config)
#' stopifnot(nrow(xx) == 80)
#' stopifnot( length(unique(xx$interaction)) == (2 + 2 + 2 * 2))
summarize_stats_factors <- function(pdata, config) {
  pdata <- complete_cases(pdata, config)
  fac_res <- list()
  stats <- summarize_stats(
    pdata,
    config,
    .completed = TRUE
  )
  fac_res[["interaction"]] <- stats

  if (config$factor_depth > 1) {
    # if 1 only then done
    for (factor in config$factor_keys_depth()) {
      stats <- summarize_stats(
        pdata,
        config,
        factor_key = factor,
        .completed = TRUE
      )
      fac_res[[factor]] <- stats
    }
  }
  intfact <- dplyr::bind_rows(fac_res)
  return(intfact)
}


#' Compute mean, sd, and CV for e.g. Peptides, or proteins, for all samples.
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @rdname summarize_stats
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' bb <- prolfqua::sim_lfq_data_protein_config()
#'
#' res1 <- summarize_stats_all(bb$data, bb$config)
#'
#' stopifnot((res1 |> dplyr::filter(group_ == "All") |> nrow()) == (res1 |> nrow()))
#' res2 <- prolfqua::sim_lfq_data_2factor_config()
#' resSt <- summarize_stats_all(res2$data, res2$config)
summarize_stats_all <- function(pdata, config, .completed = FALSE) {
  summarize_stats(pdata, config, factor_key = NULL, .completed = .completed)
}


#' summarize stats output (compute quantiles)
#' @param stats_res result of running `summarize_stats`
#' @param config AnalysisConfiguration
#' @param stats summarize either sd or CV
#' @param probs for which quantiles 10, 20 etc.
#' @rdname summarize_stats
#' @export
#' @keywords internal
#' @family stats
#' @examples
#' library(ggplot2)
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' config <- bb1$config
#' data <- bb1$data
#' stats_res <- summarize_stats(data, config)
#' sq <- summarize_stats_quantiles(stats_res, config)
#' sq <- summarize_stats_quantiles(stats_res, config, stats = "CV")
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' config <- bb$config
#' data <- bb$data
#' config$get_response()
#' stats_res <- summarize_stats(data, config)
#' sq <- summarize_stats_quantiles(stats_res, config)
#' sq <- summarize_stats_quantiles(stats_res, config, stats = "sd")
#' stats_res <- summarize_stats(data, config)
#' xx <- summarize_stats_quantiles(stats_res, config, probs = seq(0,1,by = 0.1))
#' ggplot2::ggplot(xx$long, aes(x = probs, y = quantiles, color = group_)) + geom_line() + geom_point()
#'
#'
summarize_stats_quantiles <- function(stats_res, config, stats = c("sd", "CV"), probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  stats <- match.arg(stats)
  q_column <- paste0(stats, "_quantiles")

  stats_res <- stats_res |> dplyr::filter(!is.na(!!sym(stats)))
  xx2 <- stats_res |>
    dplyr::group_by(!!!syms(config$factor_keys_depth())) |>
    tidyr::nest()

  sd_quantile_res2 <- xx2 |>
    dplyr::mutate(
      !!q_column := purrr::map(
        data,
        ~ tibble(probs = probs, quantiles = quantile(.[[stats]], probs, na.rm = TRUE))
      )
    ) |>
    dplyr::select(!!!syms(c(config$factor_keys_depth(), q_column))) |>
    tidyr::unnest(cols = dplyr::all_of(q_column))

  xx <- sd_quantile_res2 |> tidyr::unite("interaction", dplyr::all_of(config$factor_keys_depth()))
  wide <- xx |> tidyr::pivot_wider(names_from = "interaction", values_from = "quantiles")
  return(list(long = sd_quantile_res2, wide = wide))
}


.lfq_power_t_test_quantiles <- function(quantile_sd, delta = 1, min.n = 1.5, power = 0.8, sig.level = 0.05) {
  minsd <- power.t.test(delta = delta, n = min.n, sd = NULL, power = power, sig.level = sig.level)$sd
  quantile_sd <- quantile_sd |>
    mutate(sdtrimmed = dplyr::if_else(.data$quantiles < .env$minsd, .env$minsd, .data$quantiles))

  #, delta = delta, power = power, sig.level = sig.level
  get_sample_size <- function(sd) {
    power.t.test(delta = delta, sd = sd, power = power, sig.level = sig.level)$n
  }

  sample_sizes <- quantile_sd |>
    mutate(N_exact = purrr::map_dbl(!!sym("sdtrimmed"), get_sample_size), N = ceiling(!!sym("N_exact")))
  return(sample_sizes)
}
#' estimate sample sizes
#' @param quantile_sd output of `summarize_stats_quantiles`
#' @param delta effect size you are interested in
#' @param power of test
#' @param sig.level P-Value
#' @param min.n smallest n to determine
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
#' config <- bb1$config
#' data2 <- bb1$data
#' stats_res <- summarize_stats(data2, config)
#' xx <- summarize_stats_quantiles(stats_res, config, probs = c(0.5,0.8))
#' bbb <- lfq_power_t_test_quantiles_V2(xx$long)
#' bbb <- dplyr::bind_rows(bbb)
#' summary <- bbb |>
#'  dplyr::select( -N_exact, -quantiles, -sdtrimmed ) |>
#'  tidyr::pivot_wider(names_from = delta, values_from = N)
#'
lfq_power_t_test_quantiles_V2 <-
  function(quantile_sd, delta = c(0.59, 1, 2), power = 0.8, sig.level = 0.05, min.n = 1.5) {
    res <- vector(mode = "list", length = length(delta))
    for (i in seq_along(delta)) {
      res[[i]] <- .lfq_power_t_test_quantiles(
        quantile_sd,
        delta = delta[i],
        min.n = min.n,
        power = power,
        sig.level = sig.level
      )
      res[[i]]$delta <- delta[i]
    }
    res <- bind_rows(res)
    return(res)
  }


#' Compute theoretical sample sizes from factor level standard deviations
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param delta effect size you are interested in
#' @param power of test
#' @param sig.level P-Value
#' @param probs numeric vector of quantile probabilities
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#' config <- bb1$config
#' data2 <- bb1$data
#'
#' res <- lfq_power_t_test_quantiles(data2, config)
#' res$summary
#' stats_res <- summarize_stats(data2, config)
#' res <- lfq_power_t_test_quantiles(data2, config, delta = 2)
#' res <- lfq_power_t_test_quantiles(data2, config, delta = c(0.5,1,2))
#'
#'
lfq_power_t_test_quantiles <- function(
  pdata,
  config,
  delta = 1,
  power = 0.8,
  sig.level = 0.05,
  probs = seq(0.5, 0.9, by = 0.1)
) {
  if (!config$is_response_transformed) {
    warning("Intensities are not transformed yet.")
  }

  stats_res <- summarize_stats(pdata, config)
  sd <- na.omit(stats_res$sd)

  if (length(sd) > 0) {
    quantiles_sd <- quantile(sd, probs)

    sample_sizes <- expand.grid(probs = probs, delta = delta)
    quantiles_sd <- quantile(sd, sample_sizes$probs)
    sample_sizes <- add_column(sample_sizes, sd = quantiles_sd, .before = 2)
    sample_sizes <- add_column(sample_sizes, quantile = names(quantiles_sd), .before = 1)

    get_sample_size <- function(sd, delta) {
      power.t.test(delta = delta, sd = sd, power = power, sig.level = sig.level)$n
    }

    sample_sizes <- sample_sizes |> mutate(N_exact = purrr::map2_dbl(sd, delta, get_sample_size))
    sample_sizes <- sample_sizes |> mutate(N = ceiling(.data$N_exact))
    sample_sizes <- sample_sizes |> mutate(FC = round(2^delta, digits = 2))

    summary <- sample_sizes |>
      dplyr::select(-dplyr::all_of(c("N_exact", "delta"))) |>
      tidyr::pivot_wider(names_from = "FC", values_from = "N", names_prefix = "FC=")
    return(list(long = sample_sizes, summary = summary))
  } else {
    message(
      "!!! ERROR !!! No standard deviation is available,
            check if model is saturated (factor level variable).
            lfq_power_t_test_quantiles.
            !!! ERROR !!!"
    )
    return(NULL)
  }
}

#' Compute theoretical sample sizes from factor level standard deviations
#' @param stats_res data.frame `summarize_stats` output
#' @param delta effect size you are interested in
#' @param power of test
#' @param sig.level P-Value
#' @param min.n smallest n to determine
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::sim_lfq_data_peptide_config()
#'
#' ldata <- LFQData$new(bb1$data, bb1$config)
#' stats_res <- summarize_stats(ldata$data, ldata$config)
#' bb <- lfq_power_t_test_proteins(stats_res)
#'
lfq_power_t_test_proteins <- function(stats_res, delta = c(0.59, 1, 2), power = 0.8, sig.level = 0.05, min.n = 1.5) {
  stats_res <- na.omit(stats_res)
  sd_delta <- tidyr::crossing(stats_res, delta = delta)

  get_sample_size <- function(sd, delta) {
    sd_threshold <- power.t.test(delta = delta, n = min.n, sd = NULL, power = power, sig.level = sig.level)$sd
    power.t.test(delta = delta, sd = max(sd_threshold, sd), power = power, sig.level = sig.level)$n
  }
  sample_sizes <- sd_delta |> dplyr::mutate(N_exact = purrr::map2_dbl(sd, delta, get_sample_size))
  sample_sizes <- sample_sizes |> dplyr::mutate(N = ceiling(.data$N_exact))
  return(sample_sizes)
}

#' plot density distribution or ecdf of sd, mean or CV
#' @param pdata data.frame with statistics
#' @param factor_key character — factor column name for colouring
#' @param stat sd, mean or CV
#' @param ggstat either density or ecdf
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
    geom_line(stat = ggstat)
  return(p)
}
#' plot density distribution or ecdf of sd, mean or cv given intensity below and above median
#' @param pdata data.frame with statistics
#' @param factor_key character — factor column name for faceting
#' @param stat sd, mean or CV
#' @param ggstat either density or ecdf
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
    geom_line(stat = ggstat) +
    facet_wrap(factor_key)
  return(p)
}

#' plot Violin plot of sd CV or mean
#'
#' @param pdata data.frame with statistics
#' @param factor_keys_depth character vector — factor columns for grouping
#' @param stat either CV, mean or sd
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
  pdata <- pdata |> tidyr::unite("groups", factor_keys_depth)
  p <- ggplot(pdata, aes(x = .data$groups, y = .data[[stat]])) +
    geom_violin() +
    ggplot2::stat_summary(fun = median, geom = "point", size = 1, color = "black")
  return(p)
}
#' plot Violin plot of sd CV or mean given intensity lower or above median
#' @param pdata data.frame with statistics
#' @param factor_key character — factor column name for x-axis
#' @param stat either CV, mean or sd
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
  median.quartile <- function(x) {
    out <- quantile(x, probs = c(0.25, 0.5, 0.75))
    names(out) <- c("ymin", "y", "ymax")
    return(out)
  }
  pdata <- pdata |> dplyr::filter(!is.na(!!sym(stat)))

  top50 <- pdata |>
    dplyr::mutate(top = ifelse(meanAbundance > median(meanAbundance, na.rm = TRUE), "top 50", "bottom 50"))

  p <- ggplot(top50, aes(x = .data[[factor_key]], y = .data[[stat]])) +
    geom_violin() +
    stat_summary(fun = median.quartile, geom = "point", shape = 3) +
    stat_summary(fun = median, geom = "point", shape = 1) +
    facet_wrap("top")
  return(p)
}

#' plot stddev vs mean to asses stability of variance
#' @param pdata data.frame with statistics
#' @param factor_keys_depth character vector — factor columns for faceting
#' @param size how many points to sample (since scatter plot to slow for all)
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
    dplyr::sample_n(size = size) |>
    dplyr::ungroup()

  p <- ggplot(pdata, aes(x = meanAbundance, y = sd)) +
    geom_point() +
    geom_smooth(method = "loess") +
    facet_wrap(factor_keys_depth, nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}
