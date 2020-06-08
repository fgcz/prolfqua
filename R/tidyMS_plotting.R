#' visualize intensity distributions
#' @export
#' @keywords internal
#' @import ggplot2
#' @family plotting
#' @examples
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' plot_intensity_distribution_violin(LFQServiceData::sample_analysis, config)
#' analysis <- transform_work_intensity(LFQServiceData::sample_analysis, config, log2)
#' plot_intensity_distribution_violin(analysis, config)
plot_intensity_distribution_violin <- function(data, config){
  p <- ggplot(data, aes_string(x = config$table$sampleName, y = config$table$getWorkIntensity() )) +
    geom_violin() +
    theme(axis.text.x = element_text(angle = ))
  if (!config$parameter$is_intensity_transformed) {
    p <- p + scale_y_continuous(trans = 'log10')
  }
  return(p)
}

#' visualize intensity distributions
#' @export
#' @keywords internal
#' @import ggplot2
#' @family plotting
#' @rdname plot_intensity_distribution_violin
#' @examples
#'
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' plot_intensity_distribution_density(LFQServiceData::sample_analysis, config)
#' analysis <- transform_work_intensity(LFQServiceData::sample_analysis, config, log2)
#' plot_intensity_distribution_density(analysis, config)
plot_intensity_distribution_density <- function(data, config){
  p <- ggplot(data, aes_string(x = config$table$getWorkIntensity(), colour = config$table$sampleName )) +
    geom_line(stat = "density")
  if (!config$parameter$is_intensity_transformed) {
    p <- p + scale_x_continuous(trans = 'log10')
  }
  return(p)
}

#' visualize correlation among samples
#' @export
#' @keywords internal
#' @family plotting
#' @rdname plot_sample_correlation
#' @examples
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' analysis <- remove_small_intensities(LFQServiceData::sample_analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' mm <- toWideConfig(analysis, config, as.matrix = TRUE)
#' plot_sample_correlation(analysis, config)
plot_sample_correlation <- function(data, config){
  matrix <- toWideConfig(data, config, as.matrix = TRUE)$data
  M <- cor(matrix, use = "pairwise.complete.obs")
  if (nrow(M) > 12) {
    corrplot::corrplot.mixed(M,upper = "ellipse",
                             lower = "pie",
                             diag = "u",
                             tl.cex = .6,
                             tl.pos = "lt",
                             tl.col = "black",
                             mar = c(2,5,5,2))
  } else{
    corrplot::corrplot.mixed(M,upper = "ellipse",
                             lower = "number",
                             lower.col = "black",
                             tl.cex = .6,
                             number.cex = .7,
                             diag = "u",
                             tl.pos = "lt",
                             tl.col = "black",
                             mar = c(2,5,5,2))

  }
  invisible(M)
}
