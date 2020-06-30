
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
                                          log_y = FALSE,
                                          show.legend = FALSE
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
    p <- p +  geom_point(aes_string(shape = isotopeLabel), show.legend = show.legend) +
      geom_line(aes_string(linetype = isotopeLabel), show.legend = show.legend)
  }else{
    formula <- sprintf("~%s", paste(factor, collapse = " + "))
    p <- ggplot(data, aes_string(x = sample, y = intensity, group = fragment,  color = peptide))
    p <- p +  geom_point(show.legend = show.legend) + geom_line(show.legend = show.legend)
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
#' @param res data.frame
#' @param proteinName title of plot
#' @param config AnalysisConfiguration
#' @param separate if heavy and light show in one plot or with separate y axis?
#' @family aggregation
#' @family plotting
#'
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
#' bb <- setup_analysis(LFQServiceData::skylineSRM_HL_data,LFQServiceData::skylineconfig_HL )
#' nest <- bb %>% group_by(LFQServiceData::skylineconfig_HL$table$hkeysDepth()) %>% nest()
#' LFQService::plot_hierarchies_line(nest$data[[1]],
#'                                   "DUM",
#'                                   LFQServiceData::skylineconfig_HL,
#'                                   separate = TRUE)
#' LFQService::plot_hierarchies_line(nest$data[[1]],
#' "DUM",
#' LFQServiceData::skylineconfig_HL,
#' separate = TRUE,
#' show.legend = TRUE)
plot_hierarchies_line <- function(res,
                                  proteinName,
                                  config,
                                  separate = FALSE,
                                  show.legend = FALSE){

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
    log_y = !config$parameter$is_intensity_transformed,
    show.legend = show.legend
  )
  return(res)
}


#' generates peptide level plots for all Proteins
#' @export
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @family aggregation
#' @family plotting
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
#' #TODO make it work for other hiearachy levels.
#' config$table$hierarchyDepth = 2
#' #res <- plot_hierarchies_line_df(resDataStart, config)
#'
plot_hierarchies_line_df <- function(pdata, config, show.legend = FALSE){
  factor_level <- config$table$factorDepth

  hierarchy_ID <- "hierarchy_ID"
  pdata <- pdata %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysDepth()), remove = FALSE)

  xnested <- pdata %>% dplyr::group_by_at(hierarchy_ID) %>% tidyr::nest()

  figs <- xnested %>%
    dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                              plot_hierarchies_line,
                              config = config, show.legend = show.legend ) )
  return(figs$plot)
}



#' add quantline to plot
#' @export
#' @family aggregation
#' @family plotting
#' @keywords internal
#' @examples
#' #todo
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



#' compute Tukeys median polish from peptide or precursor intensities
#' @family aggregation
#' @param x a matrix
#' @param name if TRUE returns the name of the summary column
#' @return data.frame with number of rows equal to number of columns of input matrix.
#' @export
#' @keywords internal
#'
#' @examples
#' library(tidyverse)
#' medpolishPly(name = TRUE)
#' gg <- matrix(runif(20),4,5)
#' rownames(gg) <- paste0("A",1:4)
#' colnames(gg) <- make.names(1:5)
#' gg
#' mx <- medpolishPly(gg)
#'
medpolishPly <- function(x, name = FALSE){
  if (name) {
    return("medpolish")
  }
  X <- medpolish(x,na.rm = TRUE, trace.iter = FALSE, maxiter = 10);
  res <- tibble("sampleName" = names(X$col) , medpolish = X$col + X$overall)
  res
}

.extractInt <- function(pdata, expression, feature, samples ){
  pdata <- pdata %>%
    dplyr::select_at( c( samples,
                         feature,
                         expression) ) %>%
    tidyr::spread(key = samples , value = expression) %>% LFQService:::.ExtractMatrix()
  return(pdata)
}

#' Extract intensity column in wide format
#' @export
#' @keywords internal
#' @examples
#' library(dplyr)
#'
#' skylineconfig <-  LFQServiceData::skylineconfig$clone(deep=TRUE)
#' skylineconfig$table$workIntensity <- "Area"
#' xnested <- LFQServiceData::sample_analysis %>%
#'  group_by_at(skylineconfig$table$hkeysDepth()) %>%
#'  tidyr::nest()
#' x <- xnested$data[[1]]
#' nn  <- x %>% dplyr::select( setdiff(skylineconfig$table$hierarchyKeys() ,  skylineconfig$table$hkeysDepth()) ) %>%
#'  distinct() %>% nrow()
#'
#' xx <- extractIntensities(x,skylineconfig)
#' stopifnot(dim(xx)==c(nn,22))
#'
#' # change hierarchyDepth ###################
#' conf <- skylineconfig$clone(deep=TRUE)
#' conf$table$hierarchyDepth = 1
#'
#' xnested <- LFQServiceData::sample_analysis %>%
#'  group_by_at(conf$table$hkeysDepth()) %>%
#'  tidyr::nest()
#' head(xnested)
#'
#' x <- xnested$data[[1]]
#' nn  <- x %>% dplyr::select( setdiff(skylineconfig$table$hierarchyKeys(),  skylineconfig$table$hkeysDepth()) ) %>%
#'  distinct() %>% nrow()
#'
#' xx <- extractIntensities(x,conf)
#' stopifnot(dim(xx)==c(nn,22))
#'
extractIntensities <- function(pdata, config ){
  table <- config$table
  .extractInt(pdata, table$getWorkIntensity(),
              setdiff(table$hierarchyKeys(), table$hkeysDepth()),
              table$sampleName)
}
#' medpolish dataframe
#' @param pdata data.frame
#' @param expression column name with intensities
#' @param feature column name e.g. peptide ids
#' @param samples column name e.g. sampleName
#' @return data.frame
#' @export
#' @family aggregation
#' @family plotting
#' @examples
#'
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' conf$table$hierarchyDepth = 1
#' xnested <- LFQServiceData::sample_analysis %>%
#'   group_by_at(conf$table$hkeysDepth()) %>% nest()
#'
#' feature <- setdiff(conf$table$hierarchyKeys(),  conf$table$hkeysDepth())
#' x <- xnested$data[[1]]
#' bb <- medpolishPlydf(x,
#'  expression = conf$table$getWorkIntensity(),
#'   feature = feature,
#'    samples = conf$table$sampleName)
#' LFQService:::.reestablish_condition(x,bb, conf)
#'
medpolishPlydf <- function(pdata, expression, feature, samples  ){
  bb <- .extractInt(pdata,
    expression = expression,
     feature = feature,
      samples = samples)
   medpolishPly(bb)

}
#' medpolishPlydf2
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @family aggregation
#' @export
#' @examples
#'
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' conf$table$hierarchyDepth = 1
#' xnested <- LFQServiceData::sample_analysis %>%
#'   group_by_at(conf$table$hkeysDepth()) %>% nest()
#'
#' feature <- setdiff(conf$table$hierarchyKeys(),  conf$table$hkeysDepth())
#' x <- xnested$data[[1]]
#' bb <- medpolishPlydf2(x,conf)
#' LFQService:::.reestablish_condition(x,bb, conf)
#'
medpolishPlydf_config <- function(pdata, config, name=FALSE){
  if (name) {
    return("medpolish")
  }

  feature <- setdiff(config$table$hierarchyKeys(),  config$table$hkeysDepth())
  res <- medpolishPlydf(pdata,
                 expression = config$table$getWorkIntensity(),
                 feature = feature,
                 samples = config$table$sampleName)
  return(res)
}


.summarizeRobust <- function(pdata, expression, feature , samples = "samples", maxIt = 20) {
  data <- pdata %>% select_at(c(samples, feature, expression)) %>% na.omit
  ##If there is only one 1 peptide for all samples return expression of that peptide
  expname <- paste0("mean.",expression)

  if (length(unique(data[[feature]])) == 1L) {
    data$lmrob <- data[[expression]]
    data$weights <- 1
    data <- rename(data,  !!expname := !!sym(expression))
    data <- data %>% select(-!!sym(feature))
    return(data)
  }

  ## modelmatrix breaks on factors with 1 level so make vector of ones (will be intercept)
  if (length(unique(data[[samples]])) == 1L) {
    data <- data %>% group_by_at(samples) %>%
      summarize(lmrob = mean(!!sym(expression)),
                !!expname := mean(!!sym(expression)), .groups = "drop")
    data$weights <- 1
    return(data)
  }

  ## sum contrast on peptide level so sample effect will be mean over all peptides instead of reference level
  formula <- as.formula(paste0("~ -1 + " , samples, " + ", feature ))
  contr.arg <- list('contr.sum')
  names(contr.arg) <- feature
  X = model.matrix(formula , data = data, contrasts.arg = contr.arg)
  ## MASS::rlm breaks on singular values.
  ## check with base lm if singular values are present.
  ## if so, these coefficients will be zero, remove this column from model matrix
  ## rinse and repeat on reduced modelmatrix untill no singular values are present
  y <- data[[expression]]
  repeat {
    fit = .lm.fit(X,y)
    id = fit$coefficients != 0
    X = X[ , id, drop = FALSE]
    if (!any(!id)) break
  }
  ## Last step is always rlm if X > has some columns left
  if (ncol(X) > 0) {
    fit = MASS::rlm(X, y, maxit = maxIt)
    data$residuals <- fit$residuals
    usamples <- unique(data[[samples]])
    coefNames <- paste0(samples, usamples)
    lmrob <- tibble(!!samples := usamples, lmrob = fit$coefficients[coefNames])

    sumdata <- data %>%
      select(-!!sym(feature)) %>%
      group_by_at(samples) %>%
      dplyr::summarize(!!expname := mean(!!sym(expression)),
                       weights = 1 / mean(residuals^2), .groups = "drop")
    if (any(is.infinite(sumdata$weights) | is.na(sumdata$weights) | sumdata$weights > 10e6)) {
      sumdata$weights = 1
    }

    res <- inner_join(sumdata, lmrob, by = samples)
    if (sum(is.na(res[[expname]])) < sum(is.na(res$lmrob))) {
      res$lmrob <- res[[expname]]
    }
  } else {
    sumdata <- data %>%
      select(-!!sym(feature)) %>%
      group_by_at(samples) %>%
      dplyr::summarize(!!expname := mean(!!sym(expression)),
                       lmrob = mean(!!sym(expression)),
                       .groups = "drop")
    sumdata$weights <- 1
    res <- sumdata
  }
  pdata <- pdata %>% dplyr::select_at(samples) %>% distinct()
  res <- left_join(pdata, res, by = samples)
  return(res)
}



#' summarize proteins using MASS:rlm
#' @keywords internal
#' @family aggregation
#' @param pdata data
#' @param expression intensities
#' @param feature e.g. peptideIDs.
#' @param sample e.g. sampleName
#' @importFrom MASS rlm
#'
#' @examples
#' xx <- data.frame(expression = rnorm(20,0,10), feature = rep(LETTERS[1:5],4), samples= rep(letters[1:4],5))
#'
#' bb <- summarizeRobust(xx , "expression", "feature", "samples", maxIt = 20)
#' bb
#'
#' xx2 <- data.frame(log2Area = rnorm(20,0,10), peptide_Id = rep(LETTERS[1:5],4), sampleName = rep(letters[1:4],5))
#' summarizeRobust(xx2, "log2Area", "peptide_Id", "sampleName")
#'
#' summarizeRobust(checksummarizationrobust87,"log2Area", "peptide_Id", "sampleName")
#' summarizeRobust(checksummarizerobust69,"log2Area", "peptide_Id", "sampleName")
#' res <- vector(100,mode = "list")
#' for (i in 1:100) {
#'   xx3 <- xx2
#'   xx3$log2Area[sample(1:20,sample(1:15,1))] <- NA
#'   res[[i]] <- list(data = xx3, summary = summarizeRobust(xx3, "log2Area", "peptide_Id", "sampleName"))
#' }
#' summarizeRobust(xx2[xx2$peptide_Id == 'A',],"log2Area", "peptide_Id", "sampleName")
#' summarizeRobust(xx2[xx2$sampleName == 'a',],"log2Area", "peptide_Id", "sampleName")
#'
#'
#' conf <- skylineconfig$clone(deep = TRUE)
#' conf$table$hierarchyDepth = 1
#' xnested <- LFQServiceData::sample_analysis %>%
#'   group_by_at(conf$table$hkeysDepth()) %>% nest()
#'
#' feature <- setdiff(conf$table$hierarchyKeys(),  conf$table$hkeysDepth())
#' x <- xnested$data[[1]]
#' bb <- summarizeRobust(x,
#'  expression = conf$table$getWorkIntensity(),
#'   feature = feature,
#'    samples = conf$table$sampleName)
#'
#' LFQService:::.reestablish_condition(x,bb, conf)
#'
summarizeRobust <- function(pdata, expression, feature , samples, maxIt = 20) {
  pdata <- unite(pdata, "feature", all_of(feature))

  res <- .summarizeRobust(pdata, expression, feature = "feature", samples, maxIt = maxIt)
  return(res)
}

#' same as summarize robust but with config
#' @export
#' @param pdata data.frame
#' @param config AnalysisConfiguraton
#' @family aggregation
#' @keywords inernal
#' @example
#' conf <- skylineconfig$clone(deep = TRUE)
#' conf$table$hierarchyDepth = 1
#' xnested <- LFQServiceData::sample_analysis %>%
#'   group_by_at(conf$table$hkeysDepth()) %>% nest()
#'
#' feature <- setdiff(conf$table$hierarchyKeys(),  conf$table$hkeysDepth())
#' x <- xnested$data[[1]]
#' bb <- summarizeRobust2(x, conf)
#'
#' LFQService:::.reestablish_condition(x,bb, conf)
summarizeRobust_config <- function(pdata, config, name= FALSE){
  if (name) {return("lmrob")}

  feature <- setdiff(config$table$hierarchyKeys(),  config$table$hkeysDepth())
  summarizeRobust(pdata, expression = config$table$getWorkIntensity(),
                  feature = feature,
                  samples = config$table$sampleName
                  , maxIt = 20)
}




#' Summarizes the intensities within hierarchy
#'
#' @param func - a function working on a matrix of intensities for each protein.
#' @return returns list with data (data.frame) and config (AnalysisConfiguration)
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' library( LFQService )
#' library(tidyverse)
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' data <- LFQServiceData::sample_analysis
#'
#' undebug( aggregate_intensity )
#' data <- LFQService::transform_work_intensity(data, config, log2)
#' bbMed <- aggregate_intensity(data, config, .func = medpolishPlydf_config)
#'
#' bbRob <- aggregate_intensity(data, config, .func = summarizeRobust_config)
#' names(bbMed$data)
#' names(bbRob$data)
#' length(bbMed$data$medpolish)
#' length(bbRob$data$lmrob)
#' plot(bbMed$data$medpolish, bbRob$data$lmrob, log="xy", pch=".")
#' abline(0,1, col=2)
#' plot(bbMed$data$medpolish[1:100], bbRob$data$lmrob[1:100])
#' abline(0,1)
#'
aggregate_intensity <- function(data, config, .func)
{
  makeName <- .func(name = TRUE)
  config <- config$clone(deep = TRUE)

  xnested <- data %>% group_by_at(config$table$hkeysDepth()) %>% nest()

  pb <- progress::progress_bar$new(total = nrow(xnested))
  message("starting aggregation")

  res <- vector( mode = "list" , length = nrow(xnested) )
  for (i in 1:nrow(xnested)) {
    pb$tick()
    aggr <- .func(xnested$data[[i]], config)
    res[[i]] <- .reestablish_condition(xnested$data[[i]], aggr , config)
  }
  xnested[[makeName]] <- res
  #xnested <- xnested %>% dplyr::mutate(!!makeName := map2(data,
  #                                   !!sym(makeName),
  #                                   function(x, y, config){pb$tick(); }, config ))

  newconfig <- make_reduced_hierarchy_config(config,
                                             workIntensity = .func(name = TRUE),
                                             hierarchy = config$table$hkeysDepth(names = FALSE))

  unnested <- xnested %>%
    dplyr::select(config$table$hkeysDepth(), makeName) %>%
    tidyr::unnest(cols = makeName) %>%
    dplyr::ungroup()

  return(list(data = unnested, config = newconfig))
}

#' Plot feature data and result of aggretation
#' @param data data.frame before aggregation
#' @param config AnalyisConfiguration
#' @param data_aggr data.frame after aggregation
#' @param config_aggr AnalysisConfiguration of aggregated data
#' @family plotting
#' @family aggregation
#' @examples
#' library( LFQService )
#' library(tidyverse)
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' data <- LFQServiceData::sample_analysis
#'
#' undebug( aggregate_intensity )
#' data <- LFQService::transform_work_intensity(data, config, log2)
#' bbMed <- aggregate_intensity(data, config, .func = medpolishPlydf_config)
#' tmpMed <- plot_aggregation(data, config, bbMed$data, bbMed$config)
#' tmpMed$plot[[1]]
#' tmpMed$plot[[2]]
#' tmpMed$plot[[3]]
#'
#' bbRob <- aggregate_intensity(data, config, .func = summarizeRobust_config)
#' tmpRob <- plot_aggregation(data, config, bbRob$data, bbRob$config)
#' tmpRob$plot[[1]]
#' tmpRob$plot[[2]]
#' tmpRob$plot[[3]]
#'
#'
plot_aggregation <- function(data, config, data_aggr, config_reduced, show.legend= FALSE ){
  hierarchy_ID <- "hierarchy_ID"
  xnested <- data %>% group_by(!!!sym(config$table$hkeysDepth())) %>% nest()
  xnested <- xnested %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysDepth()))
  xnested_aggr <- data_aggr %>% group_by(!!!syms(config_reduced$table$hkeysDepth())) %>% nest_by(.key = "other")
  xnested_aggr <- xnested_aggr %>% tidyr::unite(hierarchy_ID , !!!syms(config$table$hkeysDepth()))
  xnested_all <- inner_join(xnested, xnested_aggr , by = hierarchy_ID )


  figs <- xnested_all %>% dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                              plot_hierarchies_line,
                              config = config,
                              show.legend = show.legend ))

  figs <- figs %>%
    dplyr::mutate(plot = map2(plot, other,
                              plot_hierarchies_add_quantline,
                              config_reduced$table$getWorkIntensity(), config ))
  figs$plot[[1]]
  return(figs)
}


#' aggregates top N intensities
#'
#' run \link{rankPrecursorsByIntensity} first
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param func function to use for aggregation
#' @param N default 3 top intensities.
#' @return list with data and new reduced configuration (config)
#' @family aggregation
#' @export
#' @keywords internal
#' @examples
#'
#'
#' library(LFQService)
#' library(tidyverse)
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' res <- removeLarge_Q_Values(LFQServiceData::spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByIntensity(res,config)
#' res %>% dplyr::select(c(config$table$hierarchyKeys(),"srm_meanInt"  ,"srm_meanIntRank")) %>%
#' distinct() %>%
#' arrange(!!!syms(c(config$table$hkeysDepth(),"srm_meanIntRank")))
#'
#' mean <- function(x){mean(x, na.rm=TRUE)}
#' sum <- function(x){sum(x, na.rm = TRUE)}
#'
#' resTOPN <- aggregateTopNIntensities(res, config, func = mean_na, N=3)
#'
#' stopifnot(dim(resTOPN$data) == c(10423, 10))
#' stopifnot( names(resTOPN) %in% c("data", "config") )
#' config$table$getWorkIntensity()
#' #debug(plot_aggregation)
#' tmpRob <- plot_aggregation(res, config, resTOPN$data, resTOPN$config, show.legen=TRUE)
#' tmpRob$plot[[4]]
#'
aggregateTopNIntensities <- function(pdata , config, func, N = 3){

  xcall <- as.list( match.call() )
  newcol <- make.names(glue::glue("srm_{deparse(xcall$func)}_{xcall$N}"))
  topInt <- pdata %>%
    dplyr::filter_at( "srm_meanIntRank", any_vars(. <= N)) %>%
    dplyr::group_by_at(c( config$table$hkeysDepth(),
                          config$table$sampleName,
                          config$table$fileName,
                          config$table$isotopeLabel,
                          config$table$factorKeys()))
  sumTopInt <- topInt %>%
    dplyr::summarize( !!newcol := func(!!sym(config$table$getWorkIntensity())),
                      ident_qValue = min(!!sym(config$table$ident_qValue))) %>%
    ungroup()

  newconfig <- make_reduced_hierarchy_config(
    config,
    workIntensity = newcol,
    hierarchy = config$table$hierarchy[1:config$table$hierarchyDepth])
  return(list(data = sumTopInt, config = newconfig))
}



#' Summarizes the intensities within hierarchy
#'
#' @param func - a function working on a matrix of intensities for each protein.
#' @return retuns function object
#' @keywords internal
#'
#' @family aggregation
#' @family deprecated
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
#' x("unnest")$data %>% dplyr::select(config$table$hierarchyKeys()[1] , "medpolish")
#' config <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' config$table$hierarchyDepth <- 1
#' x <- intensity_summary_by_hkeys(data, config, func = medpolishPly)
#'
#' x("unnest")$data
#' xnested<-x()
#' dd <- x(value = "plot")
#'
#' dd$medpolishPly[[1]]
#'
#' dd$plot[[2]]
#' # example how to add peptide count information
#'
#' tmp <- summarize_hierarchy(data, config)
#' tmp <- inner_join(tmp, x("wide")$data, by = config$table$hkeysDepth())
#'
#'
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

#' median polish from normalized peptide intensities
#' @export
#' @keywords internal
#' @family aggregation
#' @family deprecated
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




