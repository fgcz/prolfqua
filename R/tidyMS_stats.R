
#' applys func - a funciton workin on matrix for each protein and returning a vector of the same length as the number of samples
#' @export
#' @keywords internal
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep = TRUE)
#'
#' res <- summarize_cv(data, config)
#' res$CV <- res$sd/res$mean
summarize_cv <- function(data, config, all = TRUE){
  intsym <- sym(config$table$getWorkIntensity())


  hierarchyFactor <- data %>%
    dplyr::group_by(!!!syms( c(config$table$hierarchyKeys(), config$table$fkeysDepth()) )) %>%
    dplyr::summarize(n = n(),
                     not_na = sum(!is.na(!!intsym)),
                     sd = sd(!!intsym, na.rm = T),
                     mean = mean(!!intsym, na.rm = T)) %>%  dplyr::ungroup()

  hierarchyFactor <- hierarchyFactor %>% dplyr::mutate_at(config$table$fkeysDepth(), funs(as.character) )

  if (all) {
    hierarchy <- data %>%
      dplyr::group_by(!!!syms( config$table$hierarchyKeys() )) %>%
      dplyr::summarize(n = n(),
                       not_na = sum(!is.na(!!intsym)),
                       sd = sd(!!intsym,na.rm = T),
                       mean = mean(!!intsym,na.rm = T))

    hierarchy <- dplyr::mutate(hierarchy, !!config$table$factorKeys()[1] := "All")
    hierarchyFactor <- dplyr::bind_rows(hierarchyFactor,hierarchy)
  }
  if (config$parameter$is_intensity_transformed == FALSE) {
    hierarchyFactor %>% dplyr::mutate(CV = sd/mean * 100) -> hierarchyFactor
  }
  return(hierarchyFactor)
}

#' summarize stats output
#' @export
#' @keywords internal
#' @examples
#' config <- skylineconfig$clone(deep = TRUE)
#' data <- sample_analysis
#' stats_res <- summarize_cv(data, config)
#' head(stats_res)
#' summarize_cv_quantiles(stats_res, config)
#' summarize_cv_quantiles(stats_res, config, stats = "CV")
#'
#'
#' data2 <- transform_work_intensity(data, config, transformation = log2)
#' stats_res <- summarize_cv(data2, config)
#' xx <- summarize_cv_quantiles(stats_res, config, probs = seq(0,1,by = 0.1))
#' ggplot(xx$long, aes(x = probs, y = quantiles, color = Time)) + geom_line() + geom_point()
#'
summarize_cv_quantiles <- function(stats_res ,config, stats = c("sd","CV"), probs = c(0.1, 0.25, 0.5, 0.75, 0.9)){
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
  quantile_sd <- quantile_sd %>% mutate(sdtrimmed = case_when(quantiles < minsd  ~ minsd, TRUE ~ quantiles))

  #, delta = delta, power = power, sig.level = sig.level
  getSampleSize <- function(sd){
    power.t.test(delta = delta, sd = sd, power = power, sig.level = sig.level)$n
  }
  #  return(getSampleSize)

  sampleSizes <- quantile_sd %>%
    mutate( N_exact = purrr::map_dbl(sdtrimmed, getSampleSize), N = ceiling(N_exact))
  return(sampleSizes)
}
#' estimate sample sizes
#' @export
#' @keywords internal
#' @examples
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' data <- LFQServiceData::sample_analysis
#' data2 <- transform_work_intensity(data, config, transformation = log2)
#' stats_res <- summarize_cv(data2, config)
#' xx <- summarize_cv_quantiles(stats_res, config, probs = c(0.5,0.8))
#' bbb <- lfq_power_t_test_quantiles_V2(xx$long)
#' bbb <- (bind_rows(bbb))
#' summary <- bbb %>% dplyr::select( -N_exact, -quantiles, -sdtrimmed ) %>% spread(delta, N, sep="=")
#' #View(summary)
lfq_power_t_test_quantiles_V2 <- function(quantile_sd,
                                          delta = c(0.59,1,2),
                                          min.n = 1.5,
                                          power = 0.8,
                                          sig.level = 0.05){

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
#' @export
#' @keywords internal
#' @examples
#'
#' data <- LFQServiceData::sample_analysis
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#'
#' data2 <- transform_work_intensity(data, config, transformation = log2)
#'
#' res <- lfq_power_t_test_quantiles(data2, config)
#' res
#' stats_res <- summarize_cv(data2, config, all = FALSE)
#' stats_res
#' res <- lfq_power_t_test_quantiles(data2, config, delta = 2)
#' res
#' res <- lfq_power_t_test_quantiles(data2, config, delta = c(0.5,1,2))
#' res
#'
lfq_power_t_test_quantiles <- function(data,
                                       config,
                                       delta = 1,
                                       power = 0.8,
                                       sig.level = 0.05,
                                       probs = seq(0.5,0.9, by = 0.1)){

  if (!config$parameter$is_intensity_transformed) {
    warning("Intensities are not transformed yet.")
  }

  stats_res <- summarize_cv(data, config, all = FALSE)
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
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- LFQServiceData::sample_analysis
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)

#' data2 <- transform_work_intensity(data, config, transformation = log2)
#' stats_res <- summarize_cv(data2, config, all = FALSE)
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


#' applys func - a funciton workin on matrix for each protein and returning a vector of the same length as the number of samples
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep = TRUE)
#' res <- summarize_cv(data, config)
#' plot_stat_density(res, config, stat = "mean")
#' plot_stat_density(res, config, stat = "sd")
#' plot_stat_density(res, config, stat = "CV")
plot_stat_density <- function(data, config, stat = c("CV","mean","sd"), ggstat = c("density", "ecdf")){
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  p <- ggplot(data, aes_string(x = stat, colour = config$table$factorKeys()[1] )) +
    geom_line(stat = ggstat)
  return(p)
}
#' plot_stat_density_median
#' @export
#' @keywords internal
#' @examples
#'
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep = TRUE)
#' res <- summarize_cv(data, config)
#' plot_stat_density_median(res, config,"CV")
#' plot_stat_density_median(res, config,"mean")
#' plot_stat_density_median(res, config,"sd")
plot_stat_density_median <- function(data, config, stat = c("CV","mean","sd"), ggstat = c("density", "ecdf")){
  stat <- match.arg(stat)
  ggstat <- match.arg(ggstat)
  data <- data %>% dplyr::filter_at(stat, all_vars(!is.na(.)))
  res <- data %>% dplyr::mutate(top = ifelse(mean > median(mean, na.rm = TRUE),"top 50","bottom 50")) -> top50
  p <- ggplot(top50, aes_string(x = stat, colour = config$table$factorKeys()[1])) +
    geom_line(stat = ggstat) + facet_wrap("top")
  return(p)
}

#' applys func - a funciton workin on matrix for each protein and returning a vector of the same length as the number of samples
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep = TRUE)
#' res <- summarize_cv(data, config)
#' plot_stat_violin(res, config, stat = "mean")
#' plot_stat_violin(res, config, stat = "sd")
#' plot_stat_violin(res, config, stat = "CV")
plot_stat_violin <- function(data, config, stat = c("CV", "mean", "sd")){
  stat <- match.arg(stat)
  p <- ggplot(data, aes_string(x = config$table$factorKeys()[1], y = stat  )) +
    geom_violin()
  return(p)
}
#' plot_stat_violin_median
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep = TRUE)
#' res <- summarize_cv(data, config)
#' plot_stat_violin_median(res, config, stat = "mean")
plot_stat_violin_median <- function(data, config , stat = c("CV", "mean", "sd")){
  median.quartile <- function(x){
    out <- quantile(x, probs = c(0.25,0.5,0.75))
    names(out) <- c("ymin","y","ymax")
    return(out)
  }
  data <- data %>% dplyr::filter_at(stat, all_vars(!is.na(.)))

  res <- data %>%
    dplyr::mutate(top = ifelse(mean > median(mean, na.rm = TRUE),"top 50","bottom 50")) ->
    top50

  p <- ggplot(top50, aes_string(x = config$table$factorKeys()[1], y = stat)) +
    geom_violin() +
    stat_summary(fun = median.quartile, geom = 'point', shape = 3) +
    stat_summary(fun = median, geom = 'point', shape = 1) +
    facet_wrap("top")
  return(p)
}

#' stddev vs mean
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep = TRUE)
#' res <- summarize_cv(data, config)
#'
#' plot_stdv_vs_mean(res, config)
#' datalog2 <- transform_work_intensity(data, config, transformation = log2)
#' statlog2 <- summarize_cv(datalog2, config)
#' plot_stdv_vs_mean(statlog2, config)
#' config$table$getWorkIntensity()
#' config$table$popWorkIntensity()
#' datasqrt <- transform_work_intensity(data, config, transformation = sqrt)
#' ressqrt <- summarize_cv(datasqrt, config)
#' plot_stdv_vs_mean(ressqrt, config)
plot_stdv_vs_mean <- function(data, config, size=200){
  summary <- data %>%
    group_by_at(config$table$fkeysDepth()) %>%
    dplyr::summarize(n = n(),.groups="drop")
  size <- min(size, min(summary$n))

  data <- data %>%
    group_by_at(config$table$fkeysDepth()) %>%
    sample_n(size = size) %>%
    ungroup()

  p <- ggplot(data, aes(x = mean, y = abs(sd))) +
    geom_point() +
    geom_smooth(method = "loess") +
    facet_wrap(config$table$fkeysDepth(), nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}