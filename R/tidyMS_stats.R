#' compute pooled variance
#' @rdname pooled_var
#' @param x data.frame
#' @return data.frame
#' @examples
#'
#' x <- data.frame(not_na =c(1,2,2), var = c(3,4,4), mean = c(3,3,3))
#' x <- data.frame(not_na =c(1,2,1,1), var = c(NA, 0.0370, NA, NA), mean = c(-1.94,-1.46,-1.87,-1.45) )
#' prolfqua:::pooled_V2(na.omit(x))
#' prolfqua:::pooled_V1(na.omit(x))
#' x <- x[1,, drop=FALSE]
#' x
#' na.omit(x)
#' prolfqua:::pooled_V2(na.omit(x))
pooled_V2 <- function(x){

  n <- x$not_na
  sample.var <- x$var
  sample.mean <- x$mean
  pool.n <- sum(n)

  pool.mean <- sum(n * sample.mean)/pool.n
  deviation <- sample.mean - pool.mean

  SS <- (n - 1) * sample.var
  pool.SS <- sum(SS) + sum(n * deviation^2)
  pool.var <- pool.SS/(pool.n - 1)
  n.groups <-  length(sample.var)
  sdT = sqrt(pool.var * 2 / (pool.n/n.groups))

  res <- data.frame(
    n.groups = n.groups,
    n = pool.n,
    df = pool.n - n.groups,
    sd = sqrt(pool.var),
    var = pool.var,
    sdT = sdT,
    mean = pool.mean
  )
  return(res)
}

#' compute pooled variance V1
#' @rdname pooled_var
#' @param x data.frame
pooled_V1 <- function(x){
  n <- x$not_na
  sample.var <- x$var
  sample.mean <- x$mean
  pool.n <- sum(n)

  n.groups <- length(sample.var)
  SS <- (n - 1) * sample.var
  pool.var <- sum(SS)/(pool.n - n.groups)

  #SS <- (n) * sample.var
  #pool.var <- sum(SS)/(pool.n)

  pool.mean <- sum(sample.mean * n)/pool.n

  sdT = sqrt(pool.var * 2 / (pool.n/n.groups))

  res <- data.frame(
    n.groups = n.groups,
    n = pool.n,
    df = pool.n - n.groups,
    sd = sqrt(pool.var),
    sdT = sdT,
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
#' x <- data.frame(not_na =c(1,2,2), var = c(3,4,4), mean = c(3,3,3))
#' x <- data.frame(not_na =c(1,2,1,1), var = c(NA, 0.0370, NA, NA), mean = c(-1.94,-1.46,-1.87,-1.45) )
#' compute_pooled(x)
#' compute_pooled(x, method = "V2")
#' #debug(compute_pooled)
#' y <- data.frame(dilution.=c("a","b","c"),
#'      n = c(4,4,4), not_na = c(0,0,1), sd =c(NA,NA,NA),
#'      var = c(NA,NA,NA),mean = c(NaN,NaN,NaN))
#' compute_pooled(y)
#' yb <- y |> dplyr::filter(not_na > 1)
compute_pooled <- function(x, method = c("V1","V2")){
  method <- match.arg(method)
  xm <- x |> dplyr::filter(.data$not_na > 0)
  meanAll <- sum(xm$mean * xm$not_na)/sum(xm$not_na)
  not_na  = sum(xm$not_na)

  func <- pooled_V1
  if (method == "V2") {
    func <- pooled_V2
  }
  x <- x |> dplyr::filter(.data$not_na > 1)

  res <- func(x)
  if (is.na(res$mean)) {
    res$mean <- meanAll
  }
  res$meanAll <- meanAll
  res$not_na <- not_na
  return(res)
}

#' pooled variance
#' @export
#' @rdname pooled_var
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' bb <- prolfqua_data('data_ionstar')$normalized()
#' config <- old2new(bb$config)
#' data <- bb$data
#'
#' res1 <- summarize_stats(data, config)
#' pv <- poolvar(res1, config)
#' stopifnot(nrow(pv) == nrow(res1)/5)
#'
poolvar <- function(res1, config,  method = c("V1","V2")){
  method <- match.arg(method)
  resp <- res1 |> nest(data = -all_of(config$table$hierarchy_keys()) )
  pooled <- vector(length = length(resp$data), mode = "list")
  for (i in seq_along(resp$data)) {
    #print(i)
    pooled[[i]] <- compute_pooled(resp$data[[i]], method = method)
  }
  pooled =  bind_rows(pooled)
  resp$data <- NULL
  resp <- bind_cols(resp, pooled)
  resp <- resp |> mutate(!!config$table$factor_keys()[1] := "pooled")
  return(resp)
}

#' Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param all also compute for all samples (default), or only of conditions (set to FALSE)
#' @export
#' @rdname summarize_stats
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' bb <- prolfqua_data('data_ionstar')$normalized()
#' config <- old2new(bb$config)
#' data <- bb$data
#'
#' res1 <- summarize_stats(data, config)
#' d <- res1 |> dplyr::filter(protein_Id == "CON__P01030~9~NA" & peptide_Id  == "AELADQAASWLTR")
#' d <- res1 |> dplyr::filter(protein_Id == "CON__Q3SZR3~50~NA" & peptide_Id  == "EHFVDLLLSK")
#' #CON__P02769~18~NA VHKECCHGDLLECADDR
#' d <- res1 |> dplyr::filter(protein_Id == "CON__P02769~18~NA" & peptide_Id  == "VHKECCHGDLLECADDR")
#'
summarize_stats <- function(pdata, config){
  pdata <- complete_cases(pdata, config)
  intsym <- sym(config$table$get_response())
  hierarchyFactor <- pdata |>
    dplyr::group_by(!!!syms( c(config$table$hierarchy_keys(), config$table$factor_keys_depth()) )) |>
    dplyr::summarize(n = dplyr::n(),
                     not_na = sum(!is.na(!!intsym)),
                     sd = stats::sd(!!intsym, na.rm = TRUE),
                     var = stats::var(!!intsym, na.rm = TRUE),
                     mean = mean(!!intsym, na.rm = TRUE),.groups = "drop_last") |>  dplyr::ungroup()

  hierarchyFactor <- hierarchyFactor |>
    dplyr::mutate(dplyr::across(config$table$factor_keys_depth(), as.character))
  if (config$table$is_response_transformed == FALSE) {
    hierarchyFactor |> dplyr::mutate(CV = sd/mean * 100) -> hierarchyFactor
  }
  return(ungroup(hierarchyFactor))
}

#' Compute mean, sd, and CV for e.g. Peptides, or proteins, for all samples.
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param all also compute for all samples (default), or only of conditions (set to FALSE)
#' @export
#' @rdname summarize_stats
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' bb <- prolfqua_data('data_ionstar')$normalized()
#' config <- old2new(bb$config)
#' data <- bb$data
#'
#' res1 <- summarize_stats_all(data, config)
#' d <- res1 |> dplyr::filter(protein_Id == "CON__P01030~9~NA" & peptide_Id  == "AELADQAASWLTR")
#' d <- res1 |> dplyr::filter(protein_Id == "CON__Q3SZR3~50~NA" & peptide_Id  == "EHFVDLLLSK")
#' #CON__P02769~18~NA VHKECCHGDLLECADDR
#' d <- res1 |> dplyr::filter(protein_Id == "CON__P02769~18~NA" & peptide_Id  == "VHKECCHGDLLECADDR")
#' res1 |> dplyr::filter(dilution. == "pooled")
#'
summarize_stats_all <- function(pdata, config){
  pdata <- complete_cases(pdata, config)
  intsym <- sym(config$table$get_response())
  hierarchy <- pdata |>
    dplyr::group_by(!!!syms( config$table$hierarchy_keys() )) |>
    dplyr::summarize(n = dplyr::n(),
                     not_na = sum(!is.na(!!intsym)),
                     sd = sd(!!intsym,na.rm = TRUE),
                     var = sd(!!intsym,na.rm = TRUE),
                     mean = mean(!!intsym,na.rm = TRUE))

  hierarchy <- dplyr::mutate(hierarchy, !!config$table$factor_keys()[1] := "All")
  hierarchyFactor <- hierarchy
  if (config$table$is_response_transformed == FALSE) {
    hierarchyFactor |> dplyr::mutate(CV = sd/mean * 100) -> hierarchyFactor
  }
  return(ungroup(hierarchyFactor))
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
#' bb1 <- prolfqua_data('data_ionstar')$filtered()
#' config <- old2new(bb1$config$clone( deep = TRUE))
#' data <- bb1$data
#' stats_res <- summarize_stats(data, config)
#' summarize_stats_quantiles(stats_res, config)
#' summarize_stats_quantiles(stats_res, config, stats = "CV")
#'stats_res
#' bb <- prolfqua_data('data_ionstar')$normalized()
#' config <- old2new(bb$config$clone(deep = TRUE))
#' data <- bb$data
#' config$table$get_response()
#' stats_res <- summarize_stats(data, config)
#' summarize_stats_quantiles(stats_res, config)
#' summarize_stats_quantiles(stats_res, config, stats = "sd")
#'
#' stats_res <- summarize_stats(data, config)
#' xx <- summarize_stats_quantiles(stats_res, config, probs = seq(0,1,by = 0.1))
#' ggplot2::ggplot(xx$long, aes(x = probs, y = quantiles, color = dilution.)) + geom_line() + geom_point()
#'
#'
summarize_stats_quantiles <- function(stats_res,
                                   config,
                                   stats = c("sd","CV"),
                                   probs = c(0.1, 0.25, 0.5, 0.75, 0.9)){
  stats <- match.arg(stats)
  toQuantiles <- function(x, probs_i = probs) {
    tibble(probs = probs, quantiles = quantile(x, probs_i , na.rm = TRUE))
  }
  q_column <- paste0(stats,"_quantiles")


  stats_res <- stats_res |> dplyr::filter(!is.na(!!sym(stats)))
  xx2 <- stats_res |>
    dplyr::group_by(!!!syms(config$table$factor_keys_depth())) |>
    tidyr::nest()


  sd_quantile_res2 <- xx2 |>
    dplyr::mutate( !!q_column := purrr::map(data, ~toQuantiles(.[[stats]]) ))  |>
    dplyr::select(!!!syms(c(config$table$factor_keys_depth(),q_column))) |>
    tidyr::unnest(cols = c(q_column))

  xx <- sd_quantile_res2 |> tidyr::unite("interaction",config$table$factor_keys_depth())
  wide <- xx |>  spread("interaction", .data$quantiles)
  return(list(long = sd_quantile_res2, wide = wide))
}


.lfq_power_t_test_quantiles <- function(quantile_sd,
                                        delta = 1,
                                        min.n = 1.5,
                                        power = 0.8,
                                        sig.level = 0.05
){
  minsd  <- power.t.test(delta = delta,
                         n = min.n,
                         sd = NULL,
                         power = power,
                         sig.level = sig.level)$sd
  quantile_sd <- quantile_sd |>
    mutate("sdtrimmed" := case_when(quantiles < minsd  ~ minsd, TRUE ~ quantiles))

  #, delta = delta, power = power, sig.level = sig.level
  getSampleSize <- function(sd){
    power.t.test(delta = delta, sd = sd, power = power, sig.level = sig.level)$n
  }
  #  return(getSampleSize)

  sampleSizes <- quantile_sd |>
    mutate( N_exact = purrr::map_dbl(!!sym("sdtrimmed"), getSampleSize), N = ceiling(!!sym("N_exact")))
  return(sampleSizes)
}
#' estimate sample sizes
#' @param quantile_sd output of `summarize_stats_quantiles`
#' @param delta effect size you are interested in
#' @param power of test
#' @param sigma.level P-Value
#' @param min.n smallest n to determine
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' #library(ggplot2)
#'
#' bb1 <- prolfqua_data('data_ionstar')$normalized()
#' config <- old2new(bb1$config$clone( deep = TRUE))
#' data2 <- bb1$data
#' stats_res <- summarize_stats(data2, config)
#' xx <- summarize_stats_quantiles(stats_res, config, probs = c(0.5,0.8))
#' bbb <- lfq_power_t_test_quantiles_V2(xx$long)
#' bbb <- dplyr::bind_rows(bbb)
#' summary <- bbb |>
#'  dplyr::select( -N_exact, -quantiles, -sdtrimmed ) |>
#'  tidyr::spread(delta, N, sep = "=")
#' summary
lfq_power_t_test_quantiles_V2 <-
  function(quantile_sd,
           delta = c(0.59,1,2),
           power = 0.8,
           sig.level = 0.05,
           min.n = 1.5){

    res <- vector(mode = "list", length = length(delta))
    for (i in seq_along(delta)) {
      #message("i", i , "delta_i", delta[i], "\n")
      res[[i]] <- .lfq_power_t_test_quantiles(quantile_sd,
                                              delta = delta[i],
                                              min.n = min.n,
                                              power = power,
                                              sig.level = sig.level)
      res[[i]]$delta = delta[i]
    }
    res <- bind_rows(res)
    return(res)
  }


#' Compute theoretical sample sizes from factor level standard deviations
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param delta effect size you are interested in
#' @param power of test
#' @param sigma.level P-Value
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua_data('data_ionstar')$normalized()
#' config <- old2new(bb1$config$clone( deep = TRUE))
#' data2 <- bb1$data
#'
#' res <- lfq_power_t_test_quantiles(data2, config)
#' res$summary
#' stats_res <- summarize_stats(data2, config)
#' unique(stats_res$dilution.)
#' res <- lfq_power_t_test_quantiles(data2, config, delta = 2)
#' res <- lfq_power_t_test_quantiles(data2, config, delta = c(0.5,1,2))
#'
#'
lfq_power_t_test_quantiles <- function(pdata,
                                       config,
                                       delta = 1,
                                       power = 0.8,
                                       sig.level = 0.05,
                                       probs = seq(0.5,0.9, by = 0.1)){

  if (!config$table$is_response_transformed) {
    warning("Intensities are not transformed yet.")
  }

  stats_res <- summarize_stats(pdata, config)
  sd <- na.omit(stats_res$sd)

  if (length(sd) > 0) {
    quantilesSD <- quantile(sd,probs)

    sampleSizes <- expand.grid(probs = probs, delta = delta)
    quantilesSD <- quantile( sd, sampleSizes$probs )
    sampleSizes <- add_column( sampleSizes, sd = quantilesSD, .before = 2 )
    sampleSizes <- add_column( sampleSizes, quantile = names(quantilesSD), .before = 1 )

    getSampleSize <- function(sd, delta){power.t.test(delta = delta, sd = sd, power = power, sig.level = sig.level)$n}

    sampleSizes <- sampleSizes |> mutate( N_exact = purrr::map2_dbl(sd, delta, getSampleSize))
    sampleSizes <- sampleSizes |> mutate( N = ceiling(.data$N_exact))
    sampleSizes <- sampleSizes |> mutate( FC = round(2^delta, digits = 2))

    summary <- sampleSizes |> dplyr::select( -.data$N_exact, -.data$delta) |> spread(.data$FC, .data$N, sep = "=")
    return(list(long = sampleSizes, summary = summary))
  }else{
    message("!!! ERROR !!! No standard deviation is available,
            check if model is saturated (factor level variable).
            lfq_power_t_test_quantiles.
            !!! ERROR !!!")
    return(NULL)
  }
}

#' Compute theoretical sample sizes from factor level standard deviations
#' @param stats_res data.frame `summarize_stats` output
#' @param delta effect size you are interested in
#' @param power of test
#' @param sigma.level P-Value
#' @param min.n smallest n to determine
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#' bb1 <- prolfqua::prolfqua_data('data_IonstarProtein_subsetNorm')
#'
#'
#' ldata <- LFQData$new(bb1$data, old2new(bb1$config))
#' ldata <- ldata$get_sample(20)
#' stats_res <- summarize_stats(ldata$data, ldata$config)
#'
#' bb <- lfq_power_t_test_proteins(stats_res)
#'
lfq_power_t_test_proteins <- function(stats_res,
                                      delta = c(0.59,1,2),
                                      power = 0.8,
                                      sig.level = 0.05,
                                      min.n = 1.5){


  stats_res <- na.omit(stats_res)
  sd_delta <- purrr::map_df(delta, function(x){dplyr::mutate(stats_res, delta = x)} )

  getSampleSize <- function(sd, delta){
    sd_threshold <- power.t.test(delta = delta,
                                 n = min.n,
                                 sd = NULL,
                                 power = power,
                                 sig.level = sig.level)$sd
    power.t.test(delta = delta, sd = max(sd_threshold, sd), power = power, sig.level = sig.level)$n
  }
  sampleSizes <- sd_delta |> dplyr::mutate( N_exact = purrr::map2_dbl(sd, delta,  getSampleSize))
  sampleSizes <- sampleSizes |> dplyr::mutate( N = ceiling(.data$N_exact))
  return(sampleSizes)
}

#' plot density distribution or ecdf of sd, mean or CV
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param stat sd, mean or CV
#' @param ggstat either density of ecdf
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' bb1 <-prolfqua_data('data_ionstar')$filtered()
#' config <-  old2new(bb1$config$clone( deep = TRUE))
#' data <- bb1$data
#' res <- summarize_stats(data, config)
#' plot_stat_density(res, config, stat = "mean")
#' plot_stat_density(res, config, stat = "sd")
#' plot_stat_density(res, config, stat = "CV")
plot_stat_density <- function(pdata,
                              config,
                              stat = c("CV","mean","sd"),
                              ggstat = c("density", "ecdf")){
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  p <- ggplot(pdata, aes_string(x = stat, colour = config$table$factor_keys()[1] )) +
    geom_line(stat = ggstat)
  return(p)
}
#' plot density distribution or ecdf of sd, mean or cv given intensity below and above median
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param stat sd, mean or CV
#' @param ggstat either density of ecdf
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#'
#' bb1 <- prolfqua_data('data_ionstar')$filtered()
#' config <- old2new(bb1$config$clone( deep = TRUE))
#' data2 <- bb1$data
#' res <- summarize_stats(data2, config)
#' plot_stat_density_median(res, config,"CV")
#' plot_stat_density_median(res, config,"mean")
#' plot_stat_density_median(res, config,"sd")
plot_stat_density_median <- function(pdata, config, stat = c("CV","mean","sd"), ggstat = c("density", "ecdf")){
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  pdata <- pdata |> dplyr::filter(!is.na(!!sym(stat)))
  res <- pdata |> dplyr::mutate(top = ifelse(mean > median(mean, na.rm = TRUE),"top 50","bottom 50")) -> top50
  p <- ggplot(top50, aes_string(x = stat, colour = "top")) +
    geom_line(stat = ggstat) + facet_wrap(config$table$factor_keys()[1])
  return(p)
}

#' plot Violin plot of sd CV or mean
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param stat either CV, mean or sd
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' bb1 <- prolfqua_data('data_ionstar')$filtered()
#' config <- old2new(bb1$config$clone( deep = TRUE))
#' data <- bb1$data
#' res <- summarize_stats(data, config)
#' res <- summarize_stats(data, config)
#' plot_stat_violin(res, config, stat = "mean")
#' plot_stat_violin(res, config, stat = "sd")
#' plot_stat_violin(res, config, stat = "CV")
#'
plot_stat_violin <- function(pdata, config, stat = c("CV", "mean", "sd")){
  stat <- match.arg(stat)
  pdata <- pdata |> tidyr::unite("groups", config$table$factor_keys_depth())
  p <- ggplot(pdata, aes_string(x = "groups", y = stat  )) +
    geom_violin() + ggplot2::stat_summary(fun.y = median,
                                          geom = "point", size = 1, color = "black")

  return(p)
}
#' plot Violin plot of sd CV or mean given intensity lower or above median
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param stat either CV, mean or sd
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#' bb1 <- prolfqua_data('data_ionstar')$normalized()
#' config <- old2new(bb1$config$clone( deep = TRUE))
#' data <- bb1$data
#' res <- summarize_stats(data, config)
#' plot_stat_violin_median(res, config, stat = "mean")
plot_stat_violin_median <- function(pdata, config , stat = c("CV", "mean", "sd")){
  median.quartile <- function(x){
    out <- quantile(x, probs = c(0.25,0.5,0.75))
    names(out) <- c("ymin","y","ymax")
    return(out)
  }
  pdata <- pdata |> dplyr::filter(!is.na(!!sym(stat)))

  res <- pdata |>
    dplyr::mutate(top = ifelse(mean > median(mean, na.rm = TRUE),"top 50","bottom 50")) ->
    top50

  p <- ggplot(top50, aes_string(x = config$table$factor_keys()[1], y = stat)) +
    geom_violin() +
    stat_summary(fun = median.quartile, geom = 'point', shape = 3) +
    stat_summary(fun = median, geom = 'point', shape = 1) +
    facet_wrap("top")
  return(p)
}

#' plot stddev vs mean to asses stability of variance
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param size how many points to sample (since scatter plot to slow for all)
#'
#' @export
#' @keywords internal
#' @family stats
#' @examples
#'
#'
#'
#' bb1 <- prolfqua_data('data_ionstar')$filtered()
#' config <- old2new(bb1$config$clone( deep = TRUE))
#' data <- bb1$data
#' res <- summarize_stats(data, config)
#'
#' plot_stdv_vs_mean(res, config)
#' datalog2 <- transform_work_intensity(data, config, log2)
#' statlog2 <- summarize_stats(datalog2, config)
#' plot_stdv_vs_mean(statlog2, config)
#' config$table$get_response()
#' config$table$pop_response()
#' datasqrt <- transform_work_intensity(data, config, sqrt)
#' ressqrt <- summarize_stats(datasqrt, config)
#' plot_stdv_vs_mean(ressqrt, config)
#'
plot_stdv_vs_mean <- function(pdata, config, size=2000){
  summary <- pdata |>
    group_by_at(config$table$factor_keys_depth()) |>
    dplyr::summarize(n = n(),.groups = "drop")
  size <- min(size, min(summary$n))

  pdata <- pdata |>
    group_by_at(config$table$factor_keys_depth()) |>
    sample_n(size = size) |>
    ungroup()

  p <- ggplot(pdata, aes(x = mean, y = sd)) +
    geom_point() +
    geom_smooth(method = "loess") +
    facet_wrap(config$table$factor_keys_depth(), nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}
