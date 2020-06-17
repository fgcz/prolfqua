# Direct intensity manipulation ----

#' remove rows were qVal_individual_threshold exceeded
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return data.frame
#' @export
#' @keywords internal
#' @family filter functions
#' @examples
#' analysis <- LFQServiceData::spectronautDIAData250_analysis
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' res <- removeLarge_Q_Values(analysis, config)
removeLarge_Q_Values <- function(pdata, config){
  pdata <- pdata %>%
    dplyr::filter(!!sym(config$table$ident_qValue) < config$parameter$qVal_individual_threshold)
  return(pdata)
}

#' remove rows where intensity lower then threshold
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return data.frame
#' @export
#' @keywords internal
#' @family filter functions
#' @examples
#'
#' analysis <- LFQServiceData::spectronautDIAData250_analysis
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#'
#' config$table$getWorkIntensity()
#'
#' config2 <- config$clone(deep=TRUE)
#' res1 <- remove_small_intensities(analysis, config, threshold=1 )
#' res1000 <- remove_small_intensities(analysis, config2, threshold=1000 )
#' stopifnot(nrow(res1) >  nrow(res1000))
#'
remove_small_intensities <- function(pdata, config, threshold = 1){
  resData <- pdata %>% dplyr::filter(!!sym(config$table$getWorkIntensity()) >= threshold)
  return(resData)
}
#' Transform intensity
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param transformation function to transform intensities e.g. log2
#' @param intesityNewName column name for new intensity, default NULL - generates new name from name of transformation and old working intensity column name.
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' analysis <- LFQServiceData::spectronautDIAData250_analysis
#' x <- transform_work_intensity(analysis, config, transform = log2)
#' stopifnot("log2_FG.Quantity" %in% colnames(x))
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' analysis <- LFQServiceData::spectronautDIAData250_analysis
#' x <- transform_work_intensity(analysis, config, transform = asinh)
#' stopifnot("asinh_FG.Quantity" %in% colnames(x))
#'
transform_work_intensity <- function(pdata,
                                     config,
                                     transformation,
                                     intesityNewName = NULL,
                                     deep = FALSE){
  if (deep) {
    config <- config$clone(deep = TRUE)
  }
  x <- as.list( match.call() )
  if (is.null(intesityNewName)) {
    newcol <- paste(as.character(x$transformation), config$table$getWorkIntensity(), sep = "_")
  }else{
    newcol <- intesityNewName
  }

  pdata <- pdata %>% dplyr::mutate_at(config$table$getWorkIntensity(), .funs = funs(!!sym(newcol) := transformation(.)))
  config$table$setWorkIntensity(newcol)
  message("Column added : ", newcol)
  config$parameter$is_intensity_transformed = TRUE

  if (deep) {
    return( list(data = pdata, config = config) )
  } else {
    return( pdata)
  }
}



# Summarize Q Values ----
#'
#' Compute QValue summaries for each precursor
#' adds two columns srm_QValueMin - nth smallest qvalue for each precursor
#' srm_QValueNR - nr of precursors passing the threshold
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @return data.frame
#'
#' @export
#' @keywords internal
#' @param data data
#' @param config configuration
#' @examples
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' res <- summariseQValues(LFQServiceData::sample_analysis, config)
#' stopifnot(c("srm_QValueMin", "srm_QValueNR") %in% colnames(res))
#' head(res)
#' hist(unique(res$srm_QValueMin))
#' hist(unique(res$srm_QValueNR))
summariseQValues <- function(pdata,
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

  qValueSummaries <- pdata %>%
    dplyr::select(!!!syms(c(fileName, precursorIDs, config$table$ident_qValue))) %>%
    dplyr::group_by_at(precursorIDs) %>%
    dplyr::summarise_at(  c( QValue ), .funs = funs(!!QValueMin := nthbestQValue(.,qVal_minNumber_below_experiment_threshold ),
                                                    !!QValueNR  := npass(., qVal_experiment_threshold)
    ))
  pdata <- dplyr::inner_join(pdata, qValueSummaries, by = c(precursorIDs))
  message(glue::glue("Columns added : {QValueMin}, {QValueNR}"))
  return(pdata)
}

#' filter data by max and min Q Value threshold
#'
#' employs parameters ident_qValue, qVal_minNumber_below_experiment_threshold,
#' qVal_individual_threshold and qVal_experiment_threshold
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#'
#' summarize_hierarchy(LFQServiceData::sample_analysis, config)
#' res <- filter_byQValue(LFQServiceData::sample_analysis, config)
#' summarize_hierarchy(res, config)
#'
filter_byQValue <- function(pdata, config){
  data_NA <- removeLarge_Q_Values(pdata, config)
  data_NA <- summariseQValues(data_NA, config)
  data_NA_QVal <- data_NA %>%
    dplyr::filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qVal_experiment_threshold )   )
}


# Intensities to wide ----

.ExtractMatrix <- function(x, sep = "~lfq~"){
  idx <- sapply(x,is.numeric)
  xmat <- as.matrix(x[,idx])
  idcols <- x %>% dplyr::select(which(!idx == TRUE))
  if (ncol(idcols) > 0) {
    rownames(xmat) <- x %>% dplyr::select(which(!idx == TRUE)) %>%
      tidyr::unite(x, sep = sep) %>% dplyr::pull(x)
  }
  xmat
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
  pdata <- pdata %>%
    dplyr::select( c( table$sampleName,
                      setdiff(table$hierarchyKeys(),table$hkeysDepth()),
                      table$getWorkIntensity()) ) %>%
    tidyr::spread(table$sampleName, table$getWorkIntensity()) %>% .ExtractMatrix()
  return(pdata)
}

#' transform long to wide
#' @export
#' @keywords internal
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
#' @keywords internal
#' @examples
#' library(tidyverse)
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' res <- toWideConfig(LFQServiceData::sample_analysis, config)
#' res$data
#' res$annotation
#' res <- toWideConfig(LFQServiceData::sample_analysis, config, as.matrix = TRUE)
#' head(res$data)
#' res <- scale(res$data)
#'
toWideConfig <- function(data, config, as.matrix = FALSE, fileName = FALSE, sep="~lfq~"){
  if (fileName) {
    newcolname <- config$table$fileName
  }else{
    newcolname <- config$table$sampleName
  }

  ids <- dplyr::select_at(data,
                       c( config$table$sampleName, config$table$fileName, config$table$factorKeys())) %>%
    dplyr::distinct() %>% dplyr::arrange_at(newcolname)

  res <- toWide( data, c(config$table$hierarchyKeys()) ,
                 newcolname,
                 value = config$table$getWorkIntensity() )
  if (as.matrix) {
    resMat <- as.matrix(dplyr::select(res,-dplyr::one_of(config$table$hierarchyKeys())))
    names <- res %>% dplyr::select_at(config$table$hierarchyKeys()) %>%
      tidyr::unite(precursor_id, !!!syms(config$table$hierarchyKeys()), sep=sep) %>% dplyr::pull()
    rownames(resMat) <- names
    res <- resMat
  }
  return(list(data = res, annotation = ids))
}

#' make it long
#' @export
#' @keywords internal
#' @examples
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' res <- toWideConfig(LFQServiceData::sample_analysis, conf, as.matrix = TRUE)
#' res <- scale(res$data)
#'
#' xx <- gatherItBack(res,"srm_intensityScaled", conf)
#' xx <- gatherItBack(res,"srm_intensityScaled", conf, LFQServiceData::sample_analysis)
#' conf$table$getWorkIntensity() == "srm_intensityScaled"
#'
gatherItBack <- function(pdata, value, config, data = NULL, sep = "~lfq~"){
  pdata <- dplyr::bind_cols(
    tibble::tibble("row.names" := rownames(pdata)),
    tibble::as_tibble(pdata)
  )
  pdata <- tidyr::gather(pdata,key = !!config$table$sampleName, value = !!value, 2:ncol(pdata))
  pdata <- tidyr::separate(pdata, "row.names",  config$table$hierarchyKeys(), sep = sep)
  if (!is.null(data)) {
    pdata <- dplyr::inner_join(data, pdata)
    config$table$setWorkIntensity(value)
  }
  return(pdata)
}

robustscale <- function(data,
                        dim = 2,
                        center = TRUE,
                        scale = TRUE,
                        preserveScale = TRUE)
{
  medians = NULL
  if (center) {
    medians <- apply(data, dim, median, na.rm = TRUE)
    data = sweep(data, dim, medians, "-")
  }
  mads = NULL
  if (scale) {
    mads <- apply(data, dim, mad, na.rm = TRUE)
    if (preserveScale) {
      mads <- mads/mean(mads)
    }
    data = (sweep(data, dim, mads, "/"))
  }
  return(list(data = data, medians = medians, mads = mads))
}


# Functions working on Matrices go Here ----
#' robust scale warpper
#' @keywords internal
#' @export
robust_scale <- function(data){
  return(robustscale(data)$data)
}


#' apply Function To matrix
#' @export
#' @keywords internal
#' @examples
#'
#' library(tidyverse)
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' conf$table$workIntensity <- "Area"
#' res <- applyToIntensityMatrix(LFQServiceData::sample_analysis, conf, .func = base::scale)
#'
#' stopifnot("Area_base..scale" %in% colnames(res))
#' stopifnot("Area_base..scale" == conf$table$getWorkIntensity())
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' conf$table$workIntensity <- "Area"
#' res <- applyToIntensityMatrix(LFQServiceData::sample_analysis, conf$clone(deep=TRUE), .func = robust_scale)
#'
applyToIntensityMatrix <- function(data, config, .func){
  xcall <- as.list( match.call() )
  colname <- make.names( paste( config$table$getWorkIntensity(), deparse(xcall$.func), sep = "_"))
  mat <- toWideConfig(data, config, as.matrix = TRUE)$data
  mat <- .func(mat)
  data <- gatherItBack(mat, colname, config, data)
  return(data)
}

#' scale_with_subset
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' conf <- LFQServiceData::skylineconfig$clone(deep = TRUE)
#' conf$table$workIntensity <- "Area"
#'
#' sample_analysis <- LFQServiceData::sample_analysis
#' res <- transform_work_intensity(sample_analysis, conf, log2)
#' res <- scale_with_subset(res, res, conf)
#'
#'
scale_with_subset <- function(data, subset, config){
  colname <- make.names( paste( config$table$getWorkIntensity(), "subset_scaled", sep = "_"))
  subset <- toWideConfig(subset, config, as.matrix = TRUE)$data
  scales <- LFQService:::robustscale(subset)
  mat <- toWideConfig(data, config, as.matrix = TRUE)$data
  mat = sweep(mat, 2, scales$medians, "-")
  mat = sweep(mat, 2, scales$mads, "/")

  data <- gatherItBack(mat, colname, config, data)
  return(data)
}



# uncorrelation analysis ----
.findDecorrelated <- function(res, threshold = 0.65){
  if (is.null(res))
    return(NULL)
  nrtrans <- ncol(res)
  ids <- rowSums(res < threshold, na.rm = TRUE)
  names(which((nrtrans - 1) == ids))
}

.decorelatedPly <- function(pdata, corThreshold = 0.7){
  res <- LFQService::transitionCorrelationsJack(pdata)
  decorelated <- .findDecorrelated(res,threshold = corThreshold)
  tibble( row = rownames(res), srm_decorelated = rownames(res) %in% decorelated)
}

#' marks uncorrelated elements
#' @export
#' @keywords internal
#' @section TODO: do something with warnings of type "the standard deviation is zero".
#' @section TODO: do investigate In max(x, na.rm = TRUE) : no non-missing arguments to max; returning -Inf
#' @examples
#' library(LFQService)
#' data <- LFQServiceData::sample_analysis
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' data <- complete_cases(data, config)
#' mean(is.na(data$Area))
#' #debug(markDecorrelated)
#' dataI <- markDecorrelated(data, config)
#' head(dataI)
markDecorrelated <- function(data , config, minCorrelation = 0.7){
  qvalFiltX <- data %>%  dplyr::group_by_at(config$table$hierarchyKeys()[1]) %>% tidyr::nest()
  qvalFiltX <- qvalFiltX %>%
    dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))
  HLfigs2 <- qvalFiltX %>%
    dplyr::mutate(srmDecor = map(spreadMatrix, .decorelatedPly,  minCorrelation))
  unnest_res <- HLfigs2 %>%
    dplyr::select(config$table$hierarchyKeys()[1], "srmDecor") %>% tidyr::unnest()
  unnest_res <- unnest_res %>%
    tidyr::separate(col = "row",
                    into = config$table$hierarchyKeys()[-1],
                    sep = "~lfq~")
  qvalFiltX <- dplyr::inner_join(x = data, y = unnest_res, by = c(config$table$hierarchyKeys()) )
  return(qvalFiltX)
}


# Missing Value imputation ----

simpleImpute <- function(data){
  m <- apply(data,2, mean, na.rm = TRUE )
  res <- sweep(data,2,m,"-")
  dim(data)
  dim(res)
  resMean <- apply(res, 1, mean, na.rm = TRUE)
  resid <- matrix(replicate(length(m),resMean), nrow = length(resMean))
  imp <- sweep(resid,2,m,"+")
  res <- data
  res[is.na(res)] <- imp[is.na(res)]
  return(res)
}
#' imputation based on correlation assumption
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' data <- LFQServiceData::sample_analysis
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' data <- complete_cases(data, config)
#' dim(data)
#' mean(is.na(data$Area))
#' dataI <- impute_correlationBased(data, config)
#' stopifnot(dim(dataI) == c(dim(data)+c(0,1)))
#' mean(is.na(dataI$srm_ImputedIntensity)) == 0
#'
impute_correlationBased <- function(x , config){
  x <- complete_cases(x, config)
  nestedX <- x %>%  dplyr::group_by_at(config$table$hkeysDepth()) %>% tidyr::nest()
  nestedX <- nestedX %>% dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))

  gatherItback <- function(x,config){
    x <- dplyr::bind_cols(
      row = rownames(x),
      tibble::as_tibble(x)
    )
    tidyr::gather(x,key = !!config$table$sampleName, value = "srm_ImputedIntensity", 2:ncol(x))
  }
  nestedX <- nestedX %>% dplyr::mutate(imputed = map(spreadMatrix, simpleImpute))

  nestedX <- nestedX %>% dplyr::mutate(imputed = map(imputed, gatherItback, config))
  unnest_res <- nestedX %>% dplyr::select(config$table$hkeysDepth(), "imputed") %>% tidyr::unnest(cols = c(imputed))
  unnest_res <- unnest_res %>% tidyr::separate("row",config$table$hierarchyKeys()[-1], sep = "~lfq~" )

  qvalFiltX <- dplyr::inner_join(x, unnest_res,
                          by = c(config$table$hierarchyKeys(), config$table$sampleName) )
  config$table$setWorkIntensity("srm_ImputedIntensity")
  return(qvalFiltX)
}


.make_name_AinB <- function(levelA, levelB, prefix="nr_"){
  c_name <- paste(prefix, levelB, "_IN_", levelA, sep = "")
  return(c_name)
}

.nr_B_in_A <- function(data,
                      levelA,
                      levelB,
                      merge = TRUE){
  namA <- paste(levelA, collapse = "_")
  namB <- paste(levelB, collapse = "_")
  c_name <- .make_name_AinB(namA, namB)
  tmp <- data %>%
    dplyr::select_at(c(levelA, levelB)) %>%
    dplyr::distinct() %>%
    dplyr::group_by_at(levelA) %>%
    dplyr::summarize(!!c_name := n())
  if (!merge) {
    return(tmp)
  }
  data <- dplyr::inner_join(data, tmp, by = levelA )
  message("Column added : ", c_name)
  return(list(data = data, name = c_name))
}


#' Compute nr of B per A
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @keywords internal
#' @examples
#' library(LFQService)
#' library(tidyverse)
#' config <- LFQServiceData::skylineconfig$clone(deep=TRUE)
#' data <- LFQServiceData::sample_analysis
#' hierarchy <- config$table$hierarchyKeys()
#' res <- nr_B_in_A(data, config)
#' res$data %>% dplyr::select_at(c(config$table$hkeysDepth(),  res$name)) %>% distinct() %>%
#'   dplyr::pull() %>% table()
#'
#'
#'
#' resDataStart <- LFQServiceData::skylineSRM_HL_data
#' config <-  LFQServiceData::skylineconfig_HL$clone(deep=TRUE)
#' resDataStart <- setup_analysis(resDataStart , config)
#' nr_B_in_A(resDataStart, config)
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#' config$table$hierarchyDepth <- 2
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#'
nr_B_in_A <- function(pdata, config , merge = TRUE){
  levelA <- config$table$hkeysDepth()
  levelB <- config$table$hierarchyKeys()[length(levelA) + 1]
  .nr_B_in_A(pdata, levelA, levelB , merge = merge)
}


# Summarize Intensities by Intensity or NAs ----
.rankProteinPrecursors <-
  function(data,
           config,
           column = config$table$getWorkIntensity(),
           fun = function(x){ mean(x, na.rm = TRUE)},
           summaryColumn = "srm_meanInt",
           rankColumn = "srm_meanIntRank",
           rankFunction = function(x){ min_rank(desc(x)) }
  ){
  table <- config$table

  summaryPerPrecursor <- data %>%
    dplyr::group_by(!!!syms(table$hierarchyKeys())) %>%
     dplyr::summarize(!!summaryColumn := fun(!!sym(column)))

  groupedByProtein <- summaryPerPrecursor %>%
    dplyr::arrange(!!sym( table$hierarchyKeys()[1])) %>%
    dplyr::group_by(!!sym( table$hierarchyKeys()[1]))
  rankedBySummary <- groupedByProtein %>%
    dplyr::mutate(!!rankColumn := rankFunction(!!sym(summaryColumn)))

  data <- dplyr::inner_join(data, rankedBySummary)
  return(data)
}

#' ranks precursor - peptide by intensity.
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' res <- removeLarge_Q_Values(LFQServiceData::spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByIntensity(res,config)
#' X <-res %>% dplyr::select(c(config$table$hierarchyKeys(),
#'  srm_meanInt, srm_meanIntRank)) %>% distinct()
#' X %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1], "srm_meanIntRank"  )))
rankPrecursorsByIntensity <- function(pdata, config){
  summaryColumn <- "srm_meanInt"
  rankColumn <- "srm_meanIntRank"
  pdata <- .rankProteinPrecursors(pdata, config, column = config$table$getWorkIntensity(),
                               fun = function(x){ mean(x, na.rm = TRUE)},
                               summaryColumn = summaryColumn,
                               rankColumn = rankColumn,
                               rankFunction = function(x){min_rank(desc(x))}
  )

  message("Columns added : ", summaryColumn, " ",  rankColumn)
  return(pdata)
}

#' aggregates top N intensities
#'
#' run \link{rankPrecursorsByIntensity} first
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param func function to use for aggregation
#' @param N default 3 top intensities.
#' @return list with data and new reduced configuration (config)
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
#' res %>% dplyr::select(c(config$table$hierarchyKeys(),"srm_meanInt"  ,"srm_meanIntRank")) %>% distinct() %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1],"srm_meanIntRank")))
#' mean_na <- function(x){mean(x, na.rm=TRUE)}
#'
#' res <- aggregateTopNIntensities(res, config, func = mean_na, N=3)
#'
#' stopifnot(dim(res$data) == c(10423, 10))
#' stopifnot(names(res) %in% c("data", "config"))
#'
aggregateTopNIntensities <- function(pdata , config, func, N){
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
                      ident_qValue = min(!!sym(config$table$ident_qValue)))

  newconfig <- make_reduced_hierarchy_config(config,
                                             workIntensity = newcol,
                                             hierarchy = config$table$hierarchy[1:config$table$hierarchyDepth])
  return(list(data = sumTopInt, config = newconfig))
}

# Summarise NAs on lowest hierarchy ----

#' Ranks precursors by NAs (adds new column .NARank)
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' res <- removeLarge_Q_Values(LFQServiceData::spectronautDIAData250_analysis, config)
#' res <- rankPrecursorsByNAs(res,config)
#' colnames(res)
#' x <- res %>%
#'   dplyr::select(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(T)[1], "srm_NrNotNAs") %>%
#'   distinct() %>% dplyr::summarize(sum(srm_NrNotNAs)) %>% dplyr::pull()
#' stopifnot(sum(!is.na(res[[config$table$getWorkIntensity()[1]]])) == x)
#' res %>% dplyr::select(c(config$table$hierarchyKeys(),"srm_NrNotNAs"  ,"srm_NrNotNARank")) %>% distinct() %>% arrange(!!!syms(c(config$table$hierarchyKeys()[1],"srm_NrNotNARank")))
rankPrecursorsByNAs <- function(pdata, config){
  summaryColumn <- "srm_NrNotNAs"
  rankColumn <- "srm_NrNotNARank"
  pdata <- .rankProteinPrecursors(pdata, config,
                                column = config$table$getWorkIntensity(),
                                fun = function(x){sum(!is.na(x))},
                                summaryColumn = summaryColumn,
                                rankColumn = rankColumn,
                                rankFunction = function(x){min_rank(desc(x))}
  )
  message("Columns added : ", summaryColumn, " ",  rankColumn)
  return(pdata)
}

#' removes measurments with less than percent=60 missing values in factor_level = 1
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#'
#' rm(list=ls())
#' library(LFQService)
#' library(tidyverse)
#' config <- LFQServiceData::spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- LFQServiceData::spectronautDIAData250_analysis
#' data <- removeLarge_Q_Values(data, config)
#' hierarchy_counts(data, config)
#' res <- filter_factor_levels_by_missing(data, config,percent = 60)
#' data1 <-complete_cases(data, config)
#' hierarchy_counts(res, config)
#' summarize_hierarchy(res,config) %>%
#'  dplyr::filter(!!sym(paste0(config$table$hierarchyKeys()[2],"_n")) > 1)
#'
filter_factor_levels_by_missing <- function(pdata,
                                            config,
                                            percent = 60){
  table <- config$table
  summaryColumn = "srm_NrNotNAs"
  column <- table$getWorkIntensity()

  pdata <- complete_cases( pdata , config)
  nrNA = function(x){sum(!is.na(x))}
  summaryPerPrecursor <- pdata %>%
    dplyr::group_by(!!!syms( c(table$hierarchyKeys(), table$fkeysDepth() ))) %>%
     dplyr::summarize(!!"nr" := n(), !!summaryColumn := nrNA(!!sym(column))) %>%
    dplyr::mutate(fraction = !!sym(summaryColumn)/!!sym("nr") * 100 ) %>%  dplyr::ungroup()

  summaryPerPrecursorFiltered <- summaryPerPrecursor %>% dplyr::filter(fraction > percent)
  summaryPerPrecursorFiltered <- summaryPerPrecursorFiltered %>%
    dplyr::select(c(table$hierarchyKeys())) %>% dplyr::distinct()
  stopifnot(all(colnames(summaryPerPrecursorFiltered) %in% table$hierarchyKeys()))
  res <- summaryPerPrecursorFiltered %>% left_join(pdata)
  return( dplyr::ungroup(res))
}


