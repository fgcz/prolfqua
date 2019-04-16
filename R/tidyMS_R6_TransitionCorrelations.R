# Direct intensity manipulation ----
.setLargeQValuesToNA <- function(data,
                                 QValueColumn,
                                 intensityOld,
                                 thresholdQValue = 0.05,
                                 intensityNew = "IntensitiesWithNA"){

  thresholdF <- function(x,y, threshold = 0.05){ ifelse(x < threshold, y, NA)}
  data <- data %>%
    dplyr::mutate(!!intensityNew := thresholdF(!!!syms(c(QValueColumn ,intensityOld )), threshold = thresholdQValue))
  return(data)
}

#' sets intensities to NA if qVal_individual_threshold exceeded
#' @export
#' @family filter functions
#' @examples
#' analysis <- LFQService::spectronautDIAData250_analysis
#' config <- LFQService::spectronautDIAData250_config$clone(deep=TRUE)
#' res <- removeLarge_Q_Values(analysis, config)
removeLarge_Q_Values <- function(data, config){
  data <- data %>%
    dplyr::filter(!!sym(config$table$ident_qValue) < config$parameter$qVal_individual_threshold)
  return(data)
}

#' sets intensities smaller than threshold to NA
#' @export
#' @family filter functions
#' @examples
#'
#' analysis <- LFQService::spectronautDIAData250_analysis
#' config <- LFQService::spectronautDIAData250_config$clone(deep=TRUE)
#'
#' config$table$getWorkIntensity()
#'
#' config2 <- config$clone(deep=TRUE)
#' res1 <- remove_small_intensities(analysis, config, threshold=1 )
#' res1000 <- remove_small_intensities(analysis, config2, threshold=1000 )
#' stopifnot(nrow(res1) >  nrow(res1000))
remove_small_intensities <- function(data, config, threshold = 1){
  resData <- data %>% dplyr::filter(!!sym(config$table$getWorkIntensity()) >= threshold)
  return(resData)
}
#' Transform intensity
#' @export
#' @examples
#' library(tidyverse)
#' config <- LFQService::spectronautDIAData250_config$clone(deep=TRUE)
#' analysis <- LFQService::spectronautDIAData250_analysis
#' x <- transform_work_intensity(analysis, config, transform = log2)
#' stopifnot("log2_FG.Quantity" %in% colnames(x))
#' config <- LFQService::spectronautDIAData250_config$clone(deep=TRUE)
#' analysis <- LFQService::spectronautDIAData250_analysis
#' x <- transform_work_intensity(analysis, config, transform = asinh)
#' stopifnot("asinh_FG.Quantity" %in% colnames(x))
transform_work_intensity <- function(data,
                                     config,
                                     transformation,
                                     intesityNewName = NULL){
  x <- as.list( match.call() )
  if(is.null(intesityNewName)){
    newcol <- paste(as.character(x$transformation), config$table$getWorkIntensity(), sep="_")
  }else{
    newcol <- intesityNewName
  }

  data <- data %>% dplyr::mutate_at(config$table$getWorkIntensity(), .funs = funs(!!sym(newcol) := transformation(.)))
  config$table$setWorkIntensity(newcol)
  message("Column added : ", newcol)
  config$parameter$is_intensity_transformed = TRUE

  return(data)
}

#' visualize intensity distributions
#' @export
#' @import ggplot2
#' @family plotting
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' plot_intensity_distribution_violin(sample_analysis, config)
#' analysis <- transform_work_intensity(sample_analysis, config, log2)
#' plot_intensity_distribution_violin(analysis, config)
plot_intensity_distribution_violin <- function(data, config){
  p <- ggplot(data, aes_string(x = config$table$sampleName, y = config$table$getWorkIntensity() )) +
    geom_violin() +
    theme(axis.text.x = element_text(angle=90))
  if(!config$parameter$is_intensity_transformed){
    p <- p + scale_y_continuous(trans='log10')
  }
  return(p)
}

#' visualize intensity distributions
#' @export
#' @import ggplot2
#' @family plotting
#' @rdname plot_intensity_distribution_violin
#' @examples
#'
#' config <- skylineconfig$clone(deep=TRUE)
#' plot_intensity_distribution_density(sample_analysis, config)
#' analysis <- transform_work_intensity(sample_analysis, config, log2)
#' plot_intensity_distribution_density(analysis, config)
plot_intensity_distribution_density <- function(data, config){
  p <- ggplot(data, aes_string(x = config$table$getWorkIntensity(), colour = config$table$sampleName )) +
    geom_line(stat="density")
  if(!config$parameter$is_intensity_transformed){
    p <- p + scale_x_continuous(trans='log10')
  }
  return(p)
}

#' visualize correlation among samples
#' @export
#' @family plotting
#' @rdname plot_sample_correlation
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' analysis <- remove_small_intensities(sample_analysis, config)
#' analysis <- transform_work_intensity(analysis, config, log2)
#' mm <- toWideConfig(analysis, config, as.matrix = TRUE)
#' plot_sample_correlation(analysis, config)
plot_sample_correlation <- function(data, config){
  matrix <- toWideConfig(data, config, as.matrix = TRUE)
  M <- cor(matrix, use="pairwise.complete.obs")
  if(nrow(M)>12){
    corrplot::corrplot.mixed(M,upper="ellipse",
                             lower="pie",
                             diag="u",
                             tl.cex=.6,
                             tl.pos="lt",
                             tl.col="black",
                             mar=c(2,5,5,2))
  } else{
    corrplot::corrplot.mixed(M,upper="ellipse",
                             lower="number",
                             lower.col = "black",
                             tl.cex=.6,
                             number.cex = .7,
                             diag="u",
                             tl.pos="lt",
                             tl.col="black",
                             mar=c(2,5,5,2))

  }
  invisible(M)
}


# Summarize Q Values ----
#'
#' Compute QValue summaries for each precursor
#' adds two columns srm_QValueMin - nth smallest qvalue for each precursor
#' srm_QValueNR - nr of precursors passing the threshold
#' @export
#' @param data data
#' @param config configuration
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' res <- summariseQValues(sample_analysis, config)
#' stopifnot(c("srm_QValueMin", "srm_QValueNR") %in% colnames(res))
#' head(res)
#' hist(unique(res$srm_QValueMin))
#' hist(unique(res$srm_QValueNR))
summariseQValues <- function(data,
                             config
){
  QValueMin <- "srm_QValueMin"
  QValueNR <- "srm_QValueNR"

  precursorIDs <- config$table$hierarchyKeys()
  fileName <- config$table$fileName
  QValue  <- config$table$ident_qValue
  qVal_minNumber_below_experiment_threshold <- config$parameter$qVal_minNumber_below_experiment_threshold
  qVal_experiment_threshold <- config$parameter$qVal_experiment_threshold

  nthbestQValue <-  function(x,qVal_minNumber_below_experiment_threshold){sort(x)[qVal_minNumber_below_experiment_threshold]}
  npass <-  function(x,thresh = qVal_experiment_threshold){sum(x < thresh)}

  qValueSummaries <- data %>%
    dplyr::select(!!!syms(c(fileName, precursorIDs, config$table$ident_qValue))) %>%
    dplyr::group_by_at(precursorIDs) %>%
    dplyr::summarise_at(  c( QValue ), .funs = funs(!!QValueMin := nthbestQValue(.,qVal_minNumber_below_experiment_threshold ),
                                                    !!QValueNR  := npass(., qVal_experiment_threshold)
    ))
  data <- dplyr::inner_join(data, qValueSummaries, by=c(precursorIDs))
  message(glue::glue("Columns added : {QValueMin}, {QValueNR}"))
  return(data)
}

#' filter data by max and min Q Value threshold
#'
#' employs parameters ident_qValue, qVal_minNumber_below_experiment_threshold,
#' qVal_individual_threshold and qVal_experiment_threshold
#' @export
#' @examples
#' library(tidyverse)
#' config <- skylineconfig$clone(deep=TRUE)
#'
#' summarizeHierarchy(sample_analysis, config)
#' res <- filter_byQValue(sample_analysis, config)
#' summarizeHierarchy(res, config)
filter_byQValue <- function(data, config){
  data_NA <- removeLarge_Q_Values(data, config)
  data_NA <- summariseQValues(data_NA, config)
  data_NA_QVal <- data_NA %>%
    dplyr::filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qVal_experiment_threshold )   )
}


# Intensities to wide ----

.ExtractMatrix <- function(x){
  idx <- sapply(x,is.numeric)
  xmat <- as.matrix(x[,idx])
  rownames(xmat) <- x %>% dplyr::select(which(!idx==TRUE)) %>%
    tidyr::unite(x, sep="~") %>% dplyr::pull(x)
  xmat
}


#' Extract intensity column in wide format using lowest hierarchy as key.
#' @export
#' @examples
#' library(dplyr)
#' xnested <- sample_analysis %>%
#'  group_by_at(skylineconfig$table$hkeysLevel(FALSE)) %>%
#'  tidyr::nest()
#' x <- xnested$data[[1]]
#' nn  <- x %>% dplyr::select( skylineconfig$table$hkeysLevel(TRUE) ) %>%
#'  distinct() %>% nrow()
#' xx <- extractIntensities(x,skylineconfig)
#' stopifnot(dim(xx)==c(nn,22))
#'
#' # change hierarchyLevel ###################
#' conf <- skylineconfig$clone(deep=TRUE)
#' conf$table$hierarchyLevel = 2
#'
#' xnested <- sample_analysis %>%
#'  group_by_at(conf$table$hkeysLevel(FALSE)) %>%
#'  tidyr::nest()
#' head(xnested)
#'
#' x <- xnested$data[[1]]
#' nn  <- x %>% dplyr::select( conf$table$hkeysLevel(TRUE) ) %>%
#'  distinct() %>% nrow()
#' xx <- extractIntensities(x,conf)
#' stopifnot(dim(xx)==c(nn,22))
#'

extractIntensities <- function(x, configuration ){
  table <- configuration$table
  x <- x %>%
    dplyr::select( c( table$sampleName,
                      table$hkeysLevel(TRUE),
                      table$getWorkIntensity()) ) %>%
    tidyr::spread(table$sampleName, table$getWorkIntensity()) %>% .ExtractMatrix()
  return(x)
}

#' transform long to wide
#' @export
toWide <- function(data,
                   rowIDs ,
                   columnLabels ,
                   value
){
  wide <- data %>%
    dplyr::select_at(c(rowIDs, columnLabels, value  ))
  wide <- wide %>%
    tidyr::spread( key= columnLabels , value =  value )
  return(wide)
}

#' transform long to wide
#' @export
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' res <- toWideConfig(sample_analysis, skylineconfig)
#' res
#' res <- toWideConfig(sample_analysis, config, as.matrix = TRUE)
#' head(res)
#' res <- scale(res)
#'
toWideConfig <- function(data, config , as.matrix = FALSE, fileName = FALSE){
  if(fileName){
    newcolname <- config$table$fileName
  }else{
    newcolname <- config$table$sampleName
  }

  res <- toWide( data, c(config$table$hierarchyKeys()) ,
                 newcolname,
                 value = config$table$getWorkIntensity() )
  if(as.matrix){
    resMat <- as.matrix(dplyr::select(res,-dplyr::one_of(config$table$hierarchyKeys())))
    names <- res %>% dplyr::select_at(config$table$hierarchyKeys()) %>%
      tidyr::unite(precursor_id, !!!syms(config$table$hierarchyKeys()), sep="~") %>% dplyr::pull()
    rownames(resMat) <- names
    res <- resMat
  }
  return(res)
}

#' make it long
#' @export
#' @examples
#' conf <- skylineconfig$clone(deep = TRUE)
#' res <- toWideConfig(sample_analysis, conf, as.matrix = TRUE)
#' res <- scale(res)
#'
#' xx <- gatherItBack(res,"srm_intensityScaled",conf)
#' xx <- gatherItBack(res,"srm_intensityScaled",conf,sample_analysis)
#' conf$table$getWorkIntensity() == "srm_intensityScaled"
gatherItBack <- function(x,value,config,data = NULL){
  x <- dplyr::bind_cols(
    tibble::tibble("row.names" := rownames(x)),
    tibble::as_tibble(x)
  )
  x <- gather(x,key= !!config$table$sampleName, value = !!value, 2:ncol(x))
  x <- tidyr::separate(x, "row.names",  config$table$hierarchyKeys(), sep="~")
  if(!is.null(data)){
    x <- inner_join(data, x)
    config$table$setWorkIntensity(value)
  }
  return(x)
}

# Functions working on Matrices go Here ----
#' robust scale warpper
#' @export
robust_scale <- function(data){
  return(quantable::robustscale(data)$data)
}


#' apply Function To matrix
#' @export
#' @examples
#' library(tidyverse)
#' conf <- skylineconfig$clone(deep = TRUE)
#'
#' res <- applyToIntensityMatrix(sample_analysis, conf, .func = base::scale)
#' stopifnot("Area_base..scale" %in% colnames(res))
#' stopifnot("Area_base..scale" == conf$table$getWorkIntensity())
#'
#' res <- applyToIntensityMatrix(res, conf, .func = robust_scale)
applyToIntensityMatrix <- function(data, config, .func){
  x <- as.list( match.call() )
  colname <- make.names(paste(config$table$getWorkIntensity(), deparse(x$.func), sep="_"))
  mat <- toWideConfig(data, config, as.matrix = TRUE)
  mat <- .func(mat)
  data <- gatherItBack(mat, colname, config, data)
  return(data)
}



# Decorrelation analysis ----
.findDecorrelated <- function(res, threshold = 0.65){
  if(is.null(res))
    return(NULL)
  nrtrans <- ncol(res)
  ids <- rowSums(res < threshold, na.rm = TRUE)
  names(which((nrtrans-1)== ids))
}

#' finds decorrelated measues
#' @export
decorelatedPly <- function(x, corThreshold = 0.7){
  res <- LFQService::transitionCorrelationsJack(x)
  decorelated <- .findDecorrelated(res,threshold = corThreshold)
  tibble( row = rownames(res), srm_decorelated = rownames(res) %in% decorelated)
}

#' marks decorrelated elements
#' @export
#' @importFrom purrr map
#' @section TODO: do something with warnings of type "the standard deviation is zero".
#' @section TODO: do investigate In max(x, na.rm = TRUE) : no non-missing arguments to max; returning -Inf
#' @examples
#'
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#' data <- completeCases(data, config)
#' mean(is.na(data$Area))
#' dataI <- markDecorrelated(data, config)
#' head(dataI)
markDecorrelated <- function(data , config, minCorrelation = 0.7){
  qvalFiltX <- data %>%  group_by_at(config$table$hierarchyKeys()[1]) %>% nest()
  qvalFiltX <- qvalFiltX %>%
    dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))
  HLfigs2 <- qvalFiltX %>%
    dplyr::mutate(srmDecor = map(spreadMatrix, decorelatedPly,  minCorrelation))
  unnest_res <- HLfigs2 %>%
    dplyr::select(config$table$hierarchyKeys()[1], "srmDecor") %>% tidyr::unnest()
  unnest_res <- unnest_res %>% tidyr::separate("row", config$table$hierarchyKeys()[-1], sep="~")
  qvalFiltX <- inner_join(data, unnest_res, by=c(config$table$hierarchyKeys(), config$table$hierarchyKeys(TRUE)[1]) )
  return(qvalFiltX)
}


# Missing Value imputation ----

simpleImpute <- function(data){
  m <-apply(data,2, mean, na.rm=TRUE )
  res <- sweep(data,2,m,"-")
  dim(data)
  dim(res)
  resMean <- apply(res, 1, mean, na.rm = TRUE)
  resid <- matrix(replicate(length(m),resMean), nrow=length(resMean))
  imp <- sweep(resid,2,m,"+")
  res <- data
  res[is.na(res)] <- imp[is.na(res)]
  return(res)
}
#' imputation based on correlation assumption
#' @export
#' @importFrom purrr map
#' @examples
#' library(LFQService)
#' data <- sample_analysis
#' config <- skylineconfig$clone(deep=TRUE)
#' data <- completeCases(data, config)
#' mean(is.na(data$Area))
#' dataI <- impute_correlationBased(data, config)
#' mean(is.na(dataI$srm_ImputedIntensity))
#'
impute_correlationBased <- function(x , config){
  x <- completeCases(x, config)
  nestedX <- x %>%  group_by_at(config$table$hierarchyKeys()[1]) %>% nest()
  nestedX <- nestedX %>% dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))

  gatherItback <- function(x,config){
    x <- dplyr::bind_cols(
      row=rownames(x),
      tibble::as_tibble(x)
    )
    tidyr::gather(x,key= !!config$table$sampleName, value = "srm_ImputedIntensity", 2:ncol(x))
  }
  nestedX <- nestedX %>% dplyr::mutate(imputed = map(spreadMatrix, simpleImpute))

  nestedX <- nestedX %>% dplyr::mutate(imputed = map(imputed, gatherItback, config))
  unnest_res <- nestedX %>% dplyr::select(config$table$hierarchyKeys()[1], "imputed") %>% tidyr::unnest()
  unnest_res <- unnest_res %>% separate("row",config$table$hierarchyKeys()[-1], sep="~" )

  qvalFiltX <- inner_join(x, unnest_res,
                          by=c(config$table$hierarchyKeys(), config$table$sampleName) )
  config$table$setWorkIntensity("srm_ImputedIntensity")
  return(qvalFiltX)
}

#' @export
make_name <- function(levelA, levelB, prefix="nr_"){
  c_name <- paste(prefix ,levelB,"_by_",levelA,sep="")
  return(c_name)
}

#' Compute nr of B per A
#' @export
#' @examples
#' config <- skylineconfig$clone(deep=TRUE)
#' data <- sample_analysis
#' hierarchy <- config$table$hierarchyKeys()
#' res <- nr_B_in_A(data, hierarchy[1], hierarchy[2])
#' res %>% dplyr::select(hierarchy[1],  nr_peptide_Id_by_protein_Id) %>%
#' distinct() %>% dplyr::pull() %>% table()
nr_B_in_A <- function(data,
                      levelA,
                      levelB, merge=TRUE){
  c_name <-make_name(levelA, levelB)
  tmp <- data %>%
    dplyr::select_at(c(levelA, levelB)) %>%
    dplyr::distinct() %>%
    dplyr::group_by_at(levelA) %>%
    dplyr::summarise(!!c_name:=n())
  if(!merge){
    return(tmp)
  }
  data <- dplyr::inner_join(data, tmp, by=levelA )
  message("Column added : ", c_name)
  return(data)
}

# Summarize Intensities by Intensity or NAs ----
rankProteinPrecursors <- function(data,
                                  config,
                                  column = config$table$getWorkIntensity(),
                                  fun = function(x){ mean(x, na.rm=TRUE)},
                                  summaryColumn = "srm_meanInt",
                                  rankColumn = "srm_meanIntRank",
                                  rankFunction = function(x){min_rank(desc(x))}
){
  table <- config$table

  summaryPerPrecursor <-data %>%
    dplyr::group_by(!!!syms(table$hierarchyKeys())) %>%
    dplyr::summarise(!!summaryColumn := fun(!!sym(column)))

  groupedByProtein <- summaryPerPrecursor %>%
    dplyr::arrange(!!sym( table$hierarchyKeys()[1])) %>%
    dplyr::group_by(!!sym( table$hierarchyKeys()[1]))
  rankedBySummary <- groupedByProtein %>%
    dplyr::mutate(!!rankColumn := rankFunction(!!sym(summaryColumn)))

  data <- inner_join(data, rankedBySummary)
  return(data)
}

#' ranks precursor - peptide by intensity.
#'
#' @section TODO
#' @export
#' @examples
#'
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' res <- removeLarge_Q_Values(spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByIntensity(res,config)
#' X <-res %>% dplyr::select(c(config$table$hierarchyKeys(), srm_meanInt, srm_meanIntRank)) %>% distinct()
#' X %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1], "srm_meanIntRank"  )))
rankPrecursorsByIntensity <- function(data, config){
  summaryColumn <- "srm_meanInt"
  rankColumn <- "srm_meanIntRank"
  data<- rankProteinPrecursors(data, config, column = config$table$getWorkIntensity(),
                               fun = function(x){ mean(x, na.rm=TRUE)},
                               summaryColumn = summaryColumn,
                               rankColumn = rankColumn,
                               rankFunction = function(x){min_rank(desc(x))}
  )

  message("Columns added : ", summaryColumn, " ",  rankColumn)
  return(data)
}

#' aggregates top N intensities
#'
#' run \link{rankPrecursorsByIntensity} first
#' @export
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' res <- removeLarge_Q_Values(spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByIntensity(res,config)
#' res %>% dplyr::select(c(config$table$hierarchyKeys(),"srm_meanInt"  ,"srm_meanIntRank")) %>% distinct() %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1],"srm_meanIntRank")))
#' mean_na <- function(x){mean(x, na.rm=TRUE)}
#' res <- aggregateTopNIntensities(res, config, func = mean_na, N=3)
#' stopifnot(names(res) %in% c("data", "newconfig"))
#'
aggregateTopNIntensities <- function(data , config, func, N, hierarchy_level = 1){
  x <- as.list( match.call() )
  newcol <- make.names(glue::glue("srm_{deparse(x$func)}_{x$N}"))
  topInt <- data %>%
    dplyr::filter_at( "srm_meanIntRank", any_vars(. <= N)) %>%
    dplyr::group_by(!!!syms(c( config$table$hierarchyKeys()[1:hierarchy_level],
                               config$table$sampleName,
                               config$table$fileName,
                               config$table$isotopeLabel,
                               config$table$factorKeys())))
  sumTopInt <- topInt %>%
    dplyr::summarize( !!newcol := func(!!sym(config$table$getWorkIntensity())),
                      ident_qValue = min(!!sym(config$table$ident_qValue)))

  newconfig <- make_reduced_hierarchy_config(config,
                                             workIntensity = newcol,
                                             hierarchy = config$table$hierarchy[1:hierarchy_level])
  return(list(data = sumTopInt, newconfig = newconfig))
}

# Summarise NAs on lowest hierarchy ----

#' Ranks precursors by NAs (adds new column .NARank)
#' @export
#' @examples
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' res <- removeLarge_Q_Values(spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByNAs(res,config)
#' colnames(res)
#' x <- res %>%
#'   dplyr::select(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(T)[1], "srm_NrNotNAs") %>%
#'   distinct() %>% summarize(sum(srm_NrNotNAs)) %>% dplyr::pull()
#' stopifnot(sum(!is.na(res[[config$table$getWorkIntensity()[1]]])) == x)
#' res %>% dplyr::select(c(config$table$hierarchyKeys(),"srm_NrNotNAs"  ,"srm_NrNotNARank")) %>% distinct() %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1],"srm_NrNotNARank")))
rankPrecursorsByNAs <- function(data, config){
  summaryColumn <- "srm_NrNotNAs"
  rankColumn <- "srm_NrNotNARank"
  data <- rankProteinPrecursors(data, config,
                                column = config$table$getWorkIntensity(),
                                fun = function(x){sum(!is.na(x))},
                                summaryColumn = summaryColumn,
                                rankColumn = rankColumn,
                                rankFunction = function(x){min_rank(desc(x))}
  )
  message("Columns added : ", summaryColumn, " ",  rankColumn)
  return(data)
}

#' removes measurments with less than percent=60 missing values in factor_level = 1
#' @export
#' @examples
#'
#' rm(list=ls())
#' library(LFQService)
#' library(tidyverse)
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' data <- removeLarge_Q_Values(data, config)
#' hierarchyCounts(data, config)
#' res <- filter_factor_levels_by_missing(data, config,percent = 60)
#' data1 <-completeCases(data, config)
#' hierarchyCounts(res, config)
#' summarizeHierarchy(res,config) %>% dplyr::filter(!!sym(paste0(config$table$hierarchyKeys()[2],"_n")) > 1)
#'
filter_factor_levels_by_missing <- function(data,
                                            config,
                                            percent = 60,
                                            factor_level = 1 ){
  table <- config$table
  summaryColumn = "srm_NrNotNAs"
  column <- config$table$getWorkIntensity()

  data <- completeCases( data , config)
  nrNA = function(x){sum(!is.na(x))}
  summaryPerPrecursor <- data %>%
    dplyr::group_by(!!!syms( c(table$hierarchyKeys(), table$factorKeys()[1:factor_level] ))) %>%
    dplyr::summarise(!!"nr" := n(), !!summaryColumn := nrNA(!!sym(column))) %>%
    dplyr::mutate(fraction = !!sym(summaryColumn)/!!sym("nr") * 100 ) %>% ungroup()

  summaryPerPrecursorFiltered <- summaryPerPrecursor %>% dplyr::filter(fraction > percent)
  summaryPerPrecursorFiltered <- summaryPerPrecursorFiltered %>% dplyr::select(c(table$hierarchyKeys())) %>% distinct()
  stopifnot(all(colnames(summaryPerPrecursorFiltered) %in% table$hierarchyKeys()))
  res <- summaryPerPrecursorFiltered %>% left_join(data)
  return(ungroup(res))
}


