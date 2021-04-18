#
#
# x <- data.frame(not_na =c(1,2,2), var = c(3,4,4), mean = c(3,3,3))
# .compute_pooled(x)
.compute_pooled <- function(x){
  x <- x %>% dplyr::filter(.data$not_na > 0)
  var <- sum((x$var * (x$not_na - 1)))/(sum(x$not_na) - nrow(x))
  res <- data.frame(
    #n = sum(x$n) - nrow(x), #
    n = sum(x$not_na) - nrow(x),

    not_na  = sum(x$not_na),
    sd = sqrt(var),
    var = var,
    mean = sum(x$mean * x$not_na)/sum(x$not_na)
  )
  return(res)
}


.poolvar <- function(res1, config){
  resp <- res1 %>% nest(data = -all_of(config$table$hierarchyKeys()) )
  pooled =  purrr::map_df(resp$data, .compute_pooled )
  resp$data <- NULL
  resp <- bind_cols(resp, pooled)
  resp <- resp %>% mutate(!!config$table$factorKeys()[1] := "pooled")
  return(resp)
}

#' Compute mean, sd, and CV for all Peptides, or proteins, for all interactions and all samples.
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param all also compute for all samples (default), or only of conditions (set to FALSE)
#' @export
#' @keywords internal
#' @examples
#'
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua::data_ionstar$normalized()
#' config <- bb$config
#' data <- bb$data
#'
#' res1 <- summarize_stats(data, config, all = FALSE)
#' d <- res1 %>% dplyr::filter(protein_Id == "CON__P01030~9~NA" & peptide_Id  == "ILSLAQDQVGGSAEK")
#' res1 %>% dplyr::filter(dilution. == "pooled")
#' res2 <- summarize_stats(data, config, all = TRUE)
#'
summarize_stats <- function(pdata, config, all = TRUE){
  pdata <- complete_cases(pdata, config)
  intsym <- sym(config$table$getWorkIntensity())
  hierarchyFactor <- pdata %>%
    dplyr::group_by(!!!syms( c(config$table$hierarchyKeys(), config$table$fkeysDepth()) )) %>%
    dplyr::summarize(n = dplyr::n(),
                     not_na = sum(!is.na(!!intsym)),
                     sd = sd(!!intsym, na.rm = T),
                     var = var(!!intsym, na.rm = T),
                     mean = mean(!!intsym, na.rm = T),.groups = "drop_last") %>%  dplyr::ungroup()

  hierarchyFactor <- hierarchyFactor %>%
    dplyr::mutate(dplyr::across(config$table$fkeysDepth(), as.character))
  hierPool <- .poolvar(hierarchyFactor, config )
  hierarchyFactor <- bind_rows(hierarchyFactor, hierPool)
  if (all) {
    hierarchy <- pdata %>%
      dplyr::group_by(!!!syms( config$table$hierarchyKeys() )) %>%
      dplyr::summarize(n = dplyr::n(),
                       not_na = sum(!is.na(!!intsym)),
                       sd = sd(!!intsym,na.rm = T),
                       var = sd(!!intsym,na.rm = T),
                       mean = mean(!!intsym,na.rm = T))

    hierarchy <- dplyr::mutate(hierarchy, !!config$table$factorKeys()[1] := "All")
    hierarchyFactor <- dplyr::bind_rows(hierarchyFactor,hierarchy)
  }
  if (config$table$is_intensity_transformed == FALSE) {
    hierarchyFactor %>% dplyr::mutate(CV = sd/mean * 100) -> hierarchyFactor
  }
  return(hierarchyFactor)
}

#' summarize stats output (compute quantiles)
#' @param stats_res result of running `summarize_stats`
#' @param config AnalysisConfiguration
#' @param stats summarize either sd or CV
#' @param probs for which quantiles 10, 20 etc.
#' @export
#' @keywords internal
#' @examples
#' library(ggplot2)
#' library(prolfqua)
#' bb1 <- prolfqua::data_ionstar$filtered()
#' config <- bb1$config$clone( deep = TRUE)
#' data <- bb1$data
#' stats_res <- summarize_stats(data, config)
#' summarize_stats_quantiles(stats_res, config)
#' summarize_stats_quantiles(stats_res, config, stats = "CV")
#'stats_res
#' bb <- prolfqua::data_ionstar$normalized()
#' config <- bb$config$clone(deep = TRUE)
#' data <- bb$data
#' config$table$getWorkIntensity()
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
    tibble(probs = probs, quantiles = quantile(x, probs_i , na.rm = T))
  }
  q_column <- paste0(stats,"_quantiles")

  xx2 <- stats_res %>%
    dplyr::group_by(!!!syms(config$table$fkeysDepth())) %>%
    tidyr::nest()


  sd_quantile_res2 <- xx2 %>%
    dplyr::mutate( !!q_column := purrr::map(data, ~toQuantiles(.[[stats]]) ))  %>%
    dplyr::select(!!!syms(c(config$table$fkeysDepth(),q_column))) %>%
    tidyr::unnest(cols = c(q_column))

  xx <- sd_quantile_res2 %>% tidyr::unite("interaction",config$table$fkeysDepth())
  wide <- xx %>%  spread("interaction", quantiles)
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
  quantile_sd <- quantile_sd %>%
    mutate("sdtrimmed" := case_when(quantiles < minsd  ~ minsd, TRUE ~ quantiles))

  #, delta = delta, power = power, sig.level = sig.level
  getSampleSize <- function(sd){
    power.t.test(delta = delta, sd = sd, power = power, sig.level = sig.level)$n
  }
  #  return(getSampleSize)

  sampleSizes <- quantile_sd %>%
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
#' @examples
#'
#' library(tidyverse)
#' library(ggplot2)
#' library(prolfqua)
#' bb1 <- prolfqua::data_ionstar$normalized()
#' config <- bb1$config$clone( deep = TRUE)
#' data2 <- bb1$data
#' stats_res <- summarize_stats(data2, config)
#' xx <- summarize_stats_quantiles(stats_res, config, probs = c(0.5,0.8))
#' bbb <- lfq_power_t_test_quantiles_V2(xx$long)
#' bbb <- dplyr::bind_rows(bbb)
#' summary <- bbb %>%
#'  dplyr::select( -N_exact, -quantiles, -sdtrimmed ) %>%
#'  spread(delta, N, sep = "=")
#' summary
lfq_power_t_test_quantiles_V2 <-
  function(quantile_sd,
           delta = c(0.59,1,2),
           power = 0.8,
           sig.level = 0.05,
           min.n = 1.5){

    res <- vector(mode = "list", length = length(delta))
    for (i in 1:length(delta)) {
      cat("i", i , "delta_i", delta[i], "\n")
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
#' @examples
#'
#' bb1 <- prolfqua::data_ionstar$normalized()
#' config <- bb1$config$clone( deep = TRUE)
#' data2 <- bb1$data
#'
#' res <- lfq_power_t_test_quantiles(data2, config)
#' res$summary
#' stats_res <- summarize_stats(data2, config, all = FALSE)
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

  if (!config$table$is_intensity_transformed) {
    warning("Intensities are not transformed yet.")
  }

  stats_res <- summarize_stats(pdata, config, all = FALSE)
  sd <- na.omit(stats_res$sd)

  if (length(sd) > 0) {
    quantilesSD <- quantile(sd,probs)

    sampleSizes <- expand.grid(probs = probs, delta = delta)
    quantilesSD <- quantile( sd, sampleSizes$probs )
    sampleSizes <- add_column( sampleSizes, sd = quantilesSD, .before = 2 )
    sampleSizes <- add_column( sampleSizes, quantile = names(quantilesSD), .before = 1 )

    getSampleSize <- function(sd, delta){power.t.test(delta = delta, sd = sd, power = power, sig.level = sig.level)$n}

    sampleSizes <- sampleSizes %>% mutate( N_exact = purrr::map2_dbl(sd, delta, getSampleSize))
    sampleSizes <- sampleSizes %>% mutate( N = ceiling(N_exact))
    sampleSizes <- sampleSizes %>% mutate( FC = round(2^delta, digits = 2))

    summary <- sampleSizes %>% dplyr::select( -N_exact, -delta) %>% spread(FC, N, sep="=")
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
#' @examples
#' library(prolfqua)
#' library(tidyverse)
#' bb1 <- prolfqua::data_ionstar$normalized()
#' config <- bb1$config$clone( deep = TRUE)
#' data2 <- bb1$data
#' stats_res <- summarize_stats(data2, config, all = FALSE)
#' bb <- lfq_power_t_test_proteins(stats_res)
#' head(bb)
lfq_power_t_test_proteins <- function(stats_res,
                                      delta = c(0.59,1,2),
                                      power = 0.8,
                                      sig.level = 0.05,
                                      min.n = 1.5){


  stats_res <- na.omit(stats_res)
  sd_delta <- purrr::map_df(delta, function(x){mutate(stats_res, delta = x)} )

  getSampleSize <- function(sd, delta){
    sd_threshold <- power.t.test(delta = delta,
                                 n = min.n,
                                 sd = NULL,
                                 power = power,
                                 sig.level = sig.level)$sd
    #cat("delta",delta," sd ", sd, " sd_threshold ", max(sd_threshold, sd) , " sig.level ", sig.level, "\n" )

    power.t.test(delta = delta, sd = max(sd_threshold, sd), power = power, sig.level = sig.level)$n
  }
  sampleSizes <- sd_delta %>% mutate( N_exact = purrr::map2_dbl(sd, delta, getSampleSize))
  sampleSizes <- sampleSizes %>% mutate( N = ceiling(N_exact))
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
#' @examples
#' library(prolfqua)
#' library(tidyverse)
#' bb1 <- prolfqua::data_ionstar$filtered()
#' config <- bb1$config$clone( deep = TRUE)
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
  p <- ggplot(pdata, aes_string(x = stat, colour = config$table$factorKeys()[1] )) +
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
#' @examples
#'
#' library(prolfqua)
#' library(tidyverse)
#' bb1 <- prolfqua::data_ionstar$filtered()
#' config <- bb1$config$clone( deep = TRUE)
#' data2 <- bb1$data
#' res <- summarize_stats(data2, config)
#' plot_stat_density_median(res, config,"CV")
#' plot_stat_density_median(res, config,"mean")
#' plot_stat_density_median(res, config,"sd")
plot_stat_density_median <- function(pdata, config, stat = c("CV","mean","sd"), ggstat = c("density", "ecdf")){
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  pdata <- pdata %>% dplyr::filter_at(stat, all_vars(!is.na(.)))
  res <- pdata %>% dplyr::mutate(top = ifelse(mean > median(mean, na.rm = TRUE),"top 50","bottom 50")) -> top50
  p <- ggplot(top50, aes_string(x = stat, colour = "top")) +
    geom_line(stat = ggstat) + facet_wrap(config$table$factorKeys()[1])
  return(p)
}

#' plot Violin plot of sd CV or mean
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param stat either CV, mean or sd
#' @export
#' @keywords internal
#' @examples
#' library(prolfqua)
#' library(tidyverse)
#' bb1 <- prolfqua::data_ionstar$filtered()
#' config <- bb1$config$clone( deep = TRUE)
#' data <- bb1$data
#' res <- summarize_stats(data, config)
#' res <- summarize_stats(data, config)
#' plot_stat_violin(res, config, stat = "mean")
#' plot_stat_violin(res, config, stat = "sd")
#' plot_stat_violin(res, config, stat = "CV")
plot_stat_violin <- function(pdata, config, stat = c("CV", "mean", "sd")){
  stat <- match.arg(stat)
  p <- ggplot(pdata, aes_string(x = config$table$factorKeys()[1], y = stat  )) +
    geom_violin()
  return(p)
}
#' plot Violin plot of sd CV or mean given intensity lower or above median
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param stat either CV, mean or sd
#'
#' @export
#' @keywords internal
#' @examples
#' library(prolfqua)
#' library(tidyverse)
#' bb1 <- prolfqua::data_ionstar$normalized()
#' config <- bb1$config$clone( deep = TRUE)
#' data <- bb1$data
#' res <- summarize_stats(data, config)
#' plot_stat_violin_median(res, config, stat = "mean")
plot_stat_violin_median <- function(pdata, config , stat = c("CV", "mean", "sd")){
  median.quartile <- function(x){
    out <- quantile(x, probs = c(0.25,0.5,0.75))
    names(out) <- c("ymin","y","ymax")
    return(out)
  }
  pdata <- pdata %>% dplyr::filter_at(stat, all_vars(!is.na(.)))

  res <- pdata %>%
    dplyr::mutate(top = ifelse(mean > median(mean, na.rm = TRUE),"top 50","bottom 50")) ->
    top50

  p <- ggplot(top50, aes_string(x = config$table$factorKeys()[1], y = stat)) +
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
#' @examples
#'
#' library(prolfqua)
#' library(tidyverse)
#' bb1 <- prolfqua::data_ionstar$filtered()
#' config <- bb1$config$clone( deep = TRUE)
#' data <- bb1$data
#' res <- summarize_stats(data, config)
#'
#' plot_stdv_vs_mean(res, config)
#' datalog2 <- transform_work_intensity(data, config, log2)
#' statlog2 <- summarize_stats(datalog2, config)
#' plot_stdv_vs_mean(statlog2, config)
#' config$table$getWorkIntensity()
#' config$table$popWorkIntensity()
#' datasqrt <- transform_work_intensity(data, config, sqrt)
#' ressqrt <- summarize_stats(datasqrt, config)
#' plot_stdv_vs_mean(ressqrt, config)
#'
plot_stdv_vs_mean <- function(pdata, config, size=200){
  summary <- pdata %>%
    group_by_at(config$table$fkeysDepth()) %>%
    dplyr::summarize(n = n(),.groups = "drop")
  size <- min(size, min(summary$n))

  pdata <- pdata %>%
    group_by_at(config$table$fkeysDepth()) %>%
    sample_n(size = size) %>%
    ungroup()

  p <- ggplot(pdata, aes(x = mean, y = abs(sd))) +
    geom_point() +
    geom_smooth(method = "loess") +
    facet_wrap(config$table$fkeysDepth(), nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}
