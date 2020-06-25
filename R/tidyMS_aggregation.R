
# Functions - Plotting ----
# Plot peptide and fragments
plot_hierarchies_line_default <- function(data,
                                          proteinName,
                                          sample,
                                          intensity,
                                          peptide,
                                          fragment,
                                          factor,
                                          isotopeLabel,
                                          separate = FALSE,
                                          log_y = FALSE
) {
  if (length(isotopeLabel)) {
    if (separate) {
      formula <- paste(paste( isotopeLabel, collapse = "+"), "~", paste(factor , collapse = "+"))
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group = fragment,
                                   color = peptide
      ))
    }else{
      formula <- sprintf("~%s",paste(factor, collapse = " + "))
      data <- tidyr::unite(data, "fragment_label", fragment, isotopeLabel, remove = FALSE)
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group = "fragment_label",
                                   color = peptide
      ))
    }
    p <- p +  geom_point(aes_string(shape = isotopeLabel)) +
      geom_line(aes_string(linetype = isotopeLabel))
  }else{
    formula <- sprintf("~%s", paste(factor, collapse = " + "))
    p <- ggplot(data, aes_string(x = sample, y = intensity, group = fragment,  color = peptide))
    p <- p +  geom_point() + geom_line()
  }

  #p <- ggplot(data, aes_string(x = sample, y = intensity, group = fragment,  color= peptide, linetype = isotopeLabel))
  p <- p + facet_grid(as.formula(formula), scales = "free_x"   )
  p <- p + ggtitle(proteinName) + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")
  if (log_y) {
    p <- p + scale_y_continuous(trans = 'log10')
  }
  return(p)
}

#' extracts the relevant information from the configuration to make the plot.
#' @export
#' @keywords internal
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' xnested <- LFQServiceData::sample_analysis %>%
#'  group_by_at(conf$table$hkeysDepth()) %>% tidyr::nest()
#'
#' LFQService::plot_hierarchies_line(xnested$data[[1]], xnested$protein_Id[[1]],conf )
#'
plot_hierarchies_line <- function(res,
                                  proteinName,
                                  config,
                                  separate = FALSE){

  rev_hnames <- config$table$hierarchyKeys(TRUE)
  fragment <- rev_hnames[1]
  peptide <- rev_hnames[1]

  if (length(rev_hnames) > 2) {
    peptide <- rev_hnames[2]
  }
  res <- LFQService:::plot_hierarchies_line_default(
    res,
    proteinName = proteinName,
    sample = config$table$sampleName,
    intensity = config$table$getWorkIntensity(),
    peptide = peptide,
    fragment = fragment,
    factor = config$table$fkeysDepth(),
    isotopeLabel = config$table$isotopeLabel,
    separate = separate,
    log_y = !config$parameter$is_intensity_transformed
  )
  return(res)
}


#' generates peptide level plots for all Proteins
#' @export
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @keywords internal
#' @examples
#' library(tidyverse)
#' istar <- LFQServiceData::dataIonstarNormalizedPep
#'
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 20))
#' config <-  istar$config
#'
#' config$table$is_intensity_transformed <- FALSE
#' #debug(plot_hierarchies_line_df)
#' res <- plot_hierarchies_line_df(istar$data, config)
#' res[[1]]
#'
#' config$table$is_intensity_transformed <- TRUE
#' res <- plot_hierarchies_line_df(istar$data, config)
#' res[[1]]
#'
#' istar <- LFQServiceData::dataIonstarFilteredPep
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 20))
#' config <-  istar$config
#' res <- plot_hierarchies_line_df(istar$data, config)
#' config$table$is_intensity_transformed
#' res[[1]]
#' config$table$is_intensity_transformed <- TRUE
#' res <- plot_hierarchies_line_df(istar$data, config)
#' config$table$is_intensity_transformed
#' res[[1]]
#'
#' debug(LFQService:::plot_hierarchies_line_default)
#' #TODO make it work for other hiearachy levels.
#' config$table$hierarchyDepth = 2
#' #res <- plot_hierarchies_line_df(resDataStart, config)
#'
plot_hierarchies_line_df <- function(pdata, config){
  factor_level <- config$table$factorDepth

  hierarchy_ID <- "hierarchy_ID"
  pdata <- pdata %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysDepth()), remove = FALSE)

  xnested <- pdata %>% dplyr::group_by_at(hierarchy_ID) %>% tidyr::nest()

  figs <- xnested %>%
    dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                              plot_hierarchies_line,
                              config = config ) )
  return(figs$plot)
}



#' add quantline to plot
#' @export
#' @keywords internal
#' @examples
#'
plot_hierarchies_add_quantline <- function(p, data, aes_y,  configuration){
  table <- configuration$table
  p + geom_line(data = data,
                aes_string(x = table$sampleName , y = aes_y, group = 1),
                size = 1.3,
                color = "black",
                linetype = "dashed") +
    geom_point(data = data,
               aes_string(x = table$sampleName , y = aes_y, group = 1), color = "black", shape = 10)
}


.reestablish_condition <- function(data,
                                   medpolishRes,
                                   config
){
  table <- config$table
  xx <- data %>%  dplyr::select(c(table$sampleName,
                                  table$factorKeys(),
                                  table$fileName,
                                  table$isotopeLabel)) %>% dplyr::distinct()
  res <- dplyr::inner_join(xx,medpolishRes, by = table$sampleName)
  res
}


#' median polish from normalized peptide intensities
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' data <- LFQServiceData::dataIonstarNormalizedPep
#'
#' data$data <- data$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' res <- medpolish_protein_quants(data$data,
#' data$config )
#'
#' head(res("unnest")$data)
#'
medpolish_protein_quants <- function(data, config){
  protintensity <- LFQService::intensity_summary_by_hkeys(data ,
                                                          config,
                                                          medpolishPly)
  return(protintensity)
}


#' Summarizes the intensities within hierarchy
#'
#' @param func - a function working on a matrix of intensities for each protein.
#' @return retuns function object
#' @keywords internal
#' @export
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' data <- LFQServiceData::sample_analysis
#' x <- intensity_summary_by_hkeys(data, config, func = medpolishPly)
#'
#' res <- x("unnest")
#' x("unnest")$data %>% dplyr::select(config$table$hierarchyKeys()[1] , "medpolish") %>% tidyr::unnest()
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' config$table$hierarchyDepth <- 1
#' x <- intensity_summary_by_hkeys(data, config, func = medpolishPly)
#'
#' x("unnest")$data
#' xnested<-x()
#' dd <- x(value = "plot")
#' dd$medpolishPly[[1]]
#' dd$plot[[2]]
#' # example how to add peptide count information
#'
#' tmp <- summarize_hierarchy(data, config)
#' tmp <- inner_join(tmp, x("wide")$data, by = config$table$hkeysDepth())
#' head(tmp)
intensity_summary_by_hkeys <- function(data, config, func)
{
  x <- as.list( match.call() )
  makeName <- make.names(as.character(x$func))
  config <- config$clone(deep = TRUE)

  xnested <- data %>% group_by_at(config$table$hkeysDepth()) %>% nest()

  pb <- progress::progress_bar$new(total = 3 * nrow(xnested))
  message("starting aggregation")

  xnested <- xnested %>%
    dplyr::mutate(spreadMatrix = map(data, function(x,config){pb$tick(); extractIntensities(x, config)}, config))
  xnested <- xnested %>%
    dplyr::mutate(!!makeName := map( .data$spreadMatrix , function(x){pb$tick(); func(x)}))
  xnested <- xnested %>%
    dplyr::mutate(!!makeName := map2(data, !!sym(makeName), function(x, y, config){pb$tick(); .reestablish_condition(x,y, config) }, config ))


  res_fun <- function(value = c("nested", "unnest", "wide", "plot"), DEBUG = FALSE){

    value <- match.arg(value)
    if (DEBUG) {
      return(list(config = config, value = value, xnested = xnested  ))
    }

    newconfig <- make_reduced_hierarchy_config(config,
                                               workIntensity = func(name = TRUE),
                                               hierarchy = config$table$hkeysDepth(names = FALSE))

    if (value == "nested") {
      return(list(xnested = xnested, config = newconfig))
    }else if (value == "unnest" || value == "wide") {
      unnested <- xnested %>%
        dplyr::select(config$table$hkeysDepth(), makeName) %>%
        tidyr::unnest(cols = c(medpolishPly)) %>%
        dplyr::ungroup()

      if (value == "wide") {
        wide <- LFQService::toWideConfig(unnested, newconfig)
        wide$config <- newconfig
        return(wide)
      }
      return(list(data = unnested, config = newconfig))
    }else if (value == "plot") {
      hierarchy_ID <- "hierarchy_ID"
      xnested <- xnested %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysDepth()))
      figs <- xnested %>%
        dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                                  plot_hierarchies_line, config = config ))

      figs <- figs %>%
        dplyr::mutate(plot = map2(plot, !!sym(makeName) ,
                                  plot_hierarchies_add_quantline, func(name = TRUE), config ))
      return(figs)
    }
  }
  return(res_fun)
}

#' compute tukeys median polish from peptide or precursor intensities
#' @family matrix manipulation
#' @param name if TRUE returns the name of the summary column
#' @export
#' @keywords internal
#'
#' @examples
#' library(tidyverse)
#' medpolishPly(name = TRUE)
#' gg <- matrix(runif(20),4,5)
#' rownames(gg) <- make.names(1:4)
#' colnames(gg) <- make.names(1:5)
#' mx <- medpolishPly(gg)
#'
#' # compare it with other methods of protein inference
#' dd <- tidyr::gather(tibble::as_tibble(gg))
#' #x <- robust::lmRob(value ~ key, data = dd )
#' #pred_lmRob <- c(coef(x)[1] , coef(x)[1] + coef(x)[-1])
#' xl <- lm(value ~ key , data = dd)
#' pred_lm <- c(coef(xl)[1] , coef(xl)[1] + coef(xl)[-1])
#' xr <- MASS::rlm(value ~ key , data = dd)
#' pred_rlm <- c(coef(xr)[1] , coef(xr)[1] + coef(xr)[-1])
#'
#' xx <- cbind(medpolish = mx$medpolish, pred_lm = pred_lm,pred_rlm = pred_rlm )
#' head(xx)
#' matplot(xx, type = "l")
#'
medpolishPly <- function(x, name = FALSE){
  if (name) {
    return("medpolish")
  }
  X <- medpolish(x,na.rm = TRUE, trace.iter = FALSE, maxiter = 10);
  res <- tibble("sampleName" = names(X$col) , medpolish = X$col + X$overall)
  res
}




