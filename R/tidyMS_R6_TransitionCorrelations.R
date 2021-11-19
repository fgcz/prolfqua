# Direct intensity manipulation ----

#' remove rows were qVal_individual_threshold exceeded
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return data.frame
#' @export
#' @keywords internal
#' @family filtering
#' @examples
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#'
#'
#'
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
#' @family filtering
#' @examples
#'
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#'
#' config$table$getWorkIntensity()
#'
#' res1 <- remove_small_intensities(analysis, config, threshold=1 )
#' res1000 <- remove_small_intensities(analysis, config, threshold=1000 )
#' stopifnot(nrow(res1) >  nrow(res1000))
#'
remove_small_intensities <- function(pdata, config, threshold = 1){
  resData <- pdata %>% dplyr::filter(!!sym(config$table$getWorkIntensity()) >= threshold)
  return(resData)
}
#' Transform intensity
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param .func function to transform intensities e.g. log2
#' @param .funcname generates new name from name of transformation and old working intensity column name.
#' @param intesityNewName column name for new intensity, default NULL
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#' library(tidyverse)
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#' x <- transform_work_intensity(analysis, config, .func = log2)
#'
#' stopifnot("log2_FG.Quantity" %in% colnames(x))
#' config <- dd$config_f()
#' x <- transform_work_intensity(analysis, config, .func = asinh)
#' stopifnot("asinh_FG.Quantity" %in% colnames(x))
#'
transform_work_intensity <- function(pdata,
                                     config,
                                     .func,
                                     .funcname = NULL,
                                     intesityNewName = NULL,
                                     deep = FALSE){
  if (deep) {
    config <- config$clone(deep = TRUE)
  }
  .call <- as.list( match.call() )



  if (is.null(intesityNewName)) {
    .funcname <- if (is.null(.funcname)) {deparse(.call$.func)}else{.funcname}
    newcol <- paste(.funcname, config$table$getWorkIntensity(), sep = "_")
  }else{
    newcol <- intesityNewName
  }

  pdata <- pdata %>% dplyr::mutate_at(config$table$getWorkIntensity(),
                                      .funs = funs(!!sym(newcol) := .func(.)))

  config$table$setWorkIntensity(newcol)
  message("Column added : ", newcol)
  config$table$is_intensity_transformed = TRUE

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
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#' res <- summariseQValues(analysis, config)
#' stopifnot(c("srm_QValueMin", "srm_QValueNR") %in% colnames(res))
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
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#'
#' summarize_hierarchy(analysis, config)
#' res <- filter_byQValue(analysis, config)
#' summarize_hierarchy(res, config)
#'
filter_byQValue <- function(pdata, config){
  data_NA <- removeLarge_Q_Values(pdata, config)
  data_NA <- summariseQValues(data_NA, config)
  data_NA_QVal <- data_NA %>%
    dplyr::filter_at( "srm_QValueMin" , all_vars(. < config$parameter$qVal_experiment_threshold ))
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
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#' res <- toWideConfig(analysis, config)
#' res <- toWideConfig(analysis, config, as.matrix = TRUE)
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
      tidyr::unite("precursor_id", !!!syms(config$table$hierarchyKeys()), sep = sep) %>% dplyr::pull()
    rownames(resMat) <- names
    res <- resMat
  }
  return(list(data = res, annotation = ids))
}

#' make it long
#'
#' @param pdata (matrix)
#' @param value name of column to store values in. (see `gather`)
#' @param config AnalysisConfiguration
#' @param data lfqdata
#' @param sep separater to unite the hierarchy keys.
#' @export
#'
#' @keywords internal
#' @examples
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' conf <- dd$config_f()
#' analysis <- dd$analysis(dd$data,conf)
#' res <- toWideConfig(analysis, conf, as.matrix = TRUE)
#'
#' res <- scale(res$data)
#' xx <- gatherItBack(res,"srm_intensityScaled", conf)
#' xx <- gatherItBack(res,"srm_intensityScaled", conf,analysis)
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

#' compute median and mad on matrix
#' @keywords internal
#'
.get_robscales <- function(data,
                        dim = 2)
{
  medians <- apply(data, dim, median, na.rm = TRUE)
  data = sweep(data, dim, medians, "-")
  mads <- apply(data, dim, mad, na.rm = TRUE)
  return(list( medians = medians, mads = mads ) )
}

#' compute median and standard deviation for each sample
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#' library(tidyverse)
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config$clone(deep=TRUE)
#' sample_analysis <- bb$data
#' pepIntensityNormalized <- transform_work_intensity(sample_analysis, conf, log2)
#' s1 <- get_robscales(pepIntensityNormalized, conf)
#'
#' res <- scale_with_subset(pepIntensityNormalized, pepIntensityNormalized, conf)
#' s2 <- get_robscales(res$data, conf)
#' abs(mean(s1$mads) - mean(s2$mads)) < 0.1
#'
#'
get_robscales <- function(data, config){
  data <- toWideConfig(data, config, as.matrix = TRUE)$data
  scales <- .get_robscales(data)
  return(scales)
}



# Functions working on Matrices go Here ----
#' robust scale wrapper
#' @keywords internal
#' @family preprocessing
#' @export
robust_scale <- function(data, dim = 2, preserveMean = FALSE){
  scales <- .get_robscales(data, dim = dim)
  data = sweep(data, dim, scales$medians, "-")
  mads <- scales$mads/mean(scales$mads)
  data = (sweep(data, dim, mads, "/"))

  meanmed <- mean(scales$medians)
  addmean <- if (preserveMean) {meanmed} else {0}

  return(data + addmean)
}


#' apply Function To matrix
#' @param data data.frame
#' @param config AnalysisConfiguration
#' @param .func function
#' @param .funcname name of function (used for creating new column)
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#' library(tidyverse)
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' res <- applyToIntensityMatrix(data, conf, .func = base::scale)
#' stopifnot("peptide.intensity_base..scale" %in% colnames(res))
#' stopifnot("peptide.intensity_base..scale" == conf$table$getWorkIntensity())
#' conf <- bb$config$clone(deep=TRUE)
#' conf$table$workIntensity <- "peptide.intensity"
#' res <- applyToIntensityMatrix(data, conf$clone(deep=TRUE), .func = robust_scale)
#'
#' # Normalize data using the vsn method from bioconductor
#' is_vsn <- require("vsn")
#' if(is_vsn){
#'  res <- applyToIntensityMatrix(data, conf$clone(deep=TRUE), .func = vsn::justvsn)
#' }
#'
applyToIntensityMatrix <- function(data, config, .func, .funcname = NULL){
  .call <- as.list( match.call() )
  .funcname <- if (is.null(.funcname)) { deparse(.call$.func) } else {.funcname}
  colname <- make.names( paste( config$table$getWorkIntensity(), .funcname, sep = "_"))
  mat <- toWideConfig(data, config, as.matrix = TRUE)$data
  mat <- .func(mat)
  data <- gatherItBack(mat, colname, config, data)
  return(data)
}

#' scale_with_subset
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config$clone(deep=TRUE)
#' sample_analysis <- bb$data
#' conf$table$workIntensity <- "peptide.intensity"
#'
#' res <- transform_work_intensity(sample_analysis, conf, log2)
#' s1 <- get_robscales(res, conf)
#' res <- scale_with_subset(res, res, conf)
#' s2 <- get_robscales(res$data, conf)
#' stopifnot(abs(mean(s1$mads) - mean(s2$mads)) < 1e-6)
scale_with_subset <- function(data, subset, config, preserveMean = FALSE, get_scales = TRUE){

  colname <- make.names( paste( config$table$getWorkIntensity(), "subset_scaled", sep = "_"))
  subset <- toWideConfig(subset, config, as.matrix = TRUE)$data


  scales <- .get_robscales(subset)
  mat <- toWideConfig(data, config, as.matrix = TRUE)$data
  mat = sweep(mat, 2, scales$medians, "-")
  if (!any(scales$mads == 0)) {
    mads <- scales$mads/mean(scales$mads)
    mat = sweep(mat, 2, mads, "/")
  } else {
    warning("SKIPPING scaling step in scale_with_subset function.")
  }

  meanmed <- mean(scales$medians)
  addmean <- if (preserveMean) {meanmed} else {0}
  mat <- mat + addmean
  data <- gatherItBack(mat, colname, config, data)
  if (get_scales) {
    return(list(data = data, scales = scales))
  } else {
    return(data)
  }
}

#' scale within factor levels (e.g. use for pulldown data)
#'
#' @export
#' @keywords internal
#'
#' @family preprocessing
#' @examples
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config$clone(deep=TRUE)
#' sample_analysis <- bb$data
#' conf$table$workIntensity <- "peptide.intensity"
#'
#' res <- transform_work_intensity(sample_analysis, conf, log2)
#' res <- scale_with_subset_by_factors(res, res, conf)
#'
#'
scale_with_subset_by_factors <-  function(data, subset, config, preserveMean = TRUE){
  config <- config$clone(deep = TRUE)
  dl <- group_by(data, across(config$table$fkeysDepth())) %>% nest()
  sl <- group_by(subset, across(config$table$fkeysDepth())) %>% nest()
  cf <- config$clone(deep = T)
  cf$table$factors <- NULL
  cf$table$factorDepth <- 0
  N <- length(dl$data)
  res <- vector(mode = "list", N)
  scales <- vector(mode = "list", N)

  for (i in 1:(N - 1)) {
    tmp <- scale_with_subset(dl$data[[i]], sl$data[[i]] ,
                             cf$clone(deep = TRUE) ,
                             preserveMean = TRUE,
                             get_scales = TRUE)
    res[[i]] <- tmp$data
    scales[[i]] <- tmp$scales
  }
  tmp <- scale_with_subset(dl$data[[N]],
                           sl$data[[N]],
                           cf,
                           preserveMean = TRUE,
                           get_scales = TRUE)
  res[[N]] <- tmp$data
  scales[[N]] <- tmp$scales
  #names(scales) <- dl[[1]]
  resb <- dl
  resb$data <- res
  resb <- dplyr::ungroup( unnest(resb, cols = (names(resb))) )
  config$table$setWorkIntensity(cf$table$getWorkIntensity())
  return(list(data = resb, config = config, scales = list(mads = unlist(map(scales,"mads")), medians =  unlist(map(scales,"medians")))))
}

#' normalize data by log2 and robust scaling
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return list with data.frame (data) and updated config (config)
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' istar_data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' xx <- normalize_log2_robscale(istar_data, istar$config)
#' names(xx)
#' xx$config$table$workIntensity
#'
normalize_log2_robscale <- function(pdata, config){
  pepConfig <- config$clone(deep = TRUE)
  pepIntensityNormalized <- transform_work_intensity(pdata, pepConfig, log2)
  pepConfig$table$is_intensity_transformed = TRUE

  pepIntensityNormalized <- applyToIntensityMatrix(pepIntensityNormalized,
                                                   pepConfig,
                                                   .func = robust_scale)

  pepIntensityNormalized <- pepIntensityNormalized %>%
    dplyr::rename(transformedIntensity = pepConfig$table$getWorkIntensity())
  pepConfig$table$popWorkIntensity()
  pepConfig$table$setWorkIntensity("transformedIntensity")

  return(list(data = pepIntensityNormalized, config = pepConfig))
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
  res <- prolfqua::transitionCorrelationsJack(pdata)
  decorelated <- .findDecorrelated(res,threshold = corThreshold)
  tibble( row = rownames(res), srm_decorelated = rownames(res) %in% decorelated)
}

#' marks uncorrelated peptides
#'
#' @export
#' @keywords internal
#' @section TODO: do something with warnings of type "the standard deviation is zero".
#' @section TODO: do investigate In max(x, na.rm = TRUE) : no non-missing arguments to max; returning -Inf
#' @examples
#'
#' library(prolfqua)
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' dataI <- markDecorrelated(data, conf)
#' head(dataI)
#'
markDecorrelated <- function(data , config, minCorrelation = 0.7){
  qvalFiltX <- data %>%  dplyr::group_by_at(config$table$hierarchyKeys()[1]) %>% tidyr::nest()
  qvalFiltX <- qvalFiltX %>%
    dplyr::mutate(spreadMatrix = map(data, extractIntensities, config))
  HLfigs2 <- qvalFiltX %>%
    dplyr::mutate(srmDecor = map(.data$spreadMatrix, .decorelatedPly,  minCorrelation))
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
#' library(prolfqua)
#' library(tidyverse)
#'
#' bb <- prolfqua_data('data_ionstar')$normalized()
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' mean(is.na(data$peptide.intensity))
#' dataI <- impute_correlationBased(data, config)
#' dim(dataI)
#' stopifnot(dim(dataI) == c(dim(data)+c(0,1)))
#' stopifnot(mean(is.na(dataI$srm_ImputedIntensity)) <= mean(is.na(data$peptide.intensity)))
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

  nestedX <- nestedX %>% dplyr::mutate(imputed = map(.data$spreadMatrix, simpleImpute))

  nestedX <- nestedX %>% dplyr::mutate(imputed = map(.data$imputed, gatherItback, config))
  unnest_res <- nestedX %>% dplyr::select(config$table$hkeysDepth(), "imputed") %>% tidyr::unnest(cols = .data$imputed)
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
#' library(tidyverse)
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data %>% select(-all_of("nr_peptide_Id_IN_protein_Id"))
#' hierarchy <- config$table$hierarchyKeys()
#' res <- nr_B_in_A(data, config)
#'
#' res$data %>%
#'   dplyr::select_at(c(config$table$hkeysDepth(),  res$name)) %>%
#'   distinct() %>%
#'   dplyr::pull() %>% table()
#'
#'
#' bb <- prolfqua_data('data_skylineSRM_HL_A')
#' config <- bb$config_f()
#' resDataStart <- bb$analysis(bb$data, bb$config_f())
#'
#' nr_B_in_A(resDataStart, config)
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#' config$table$hierarchyDepth <- 2
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#'
#' bb <- prolfqua_data('data_IonstarProtein_subsetNorm')
#' nr_B_in_A(bb$data, bb$config)
#' #undebug(nr_B_in_A)
nr_B_in_A <- function(pdata, config , merge = TRUE){
  levelA <- config$table$hkeysDepth()
  levelB <- config$table$hierarchyKeys()[length(levelA) + 1]
  if (is.na(levelB)) {
    warning("here is no B in A")
    return(NULL)
  }else{
    .nr_B_in_A(pdata, levelA, levelB , merge = merge)
  }
}



#' how many peptides per protein in each sample
#' @export
#' @keywords internal
#' @family summary
#' @examples
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' nr_B_in_A_per_sample(data, configur, nested =FALSE)
#' bb <- prolfqua_data('data_IonstarProtein_subsetNorm')
#' nr_B_in_A_per_sample(bb$data, bb$config, nested=FALSE)
nr_B_in_A_per_sample <- function(data, config, nested = TRUE){
  cf <- config

  levelA <- cf$table$hkeysDepth()
  levelB <- cf$table$hierarchyKeys()[length(levelA) + 1]
  if (is.na(levelB)) {
    warning("here is no B in A")
  }
  data <- prolfqua::complete_cases(data, cf)
  data <- data %>%
    dplyr::mutate(presentabsent = case_when(!is.na(!!sym(cf$table$getWorkIntensity())) ~ 1,
                                            TRUE ~ 0))
  pepStats <- data %>% group_by_at(c(cf$table$hkeysDepth(), cf$table$sampleName)) %>%
    summarize(nrPep = n(), present = sum(.data$presentabsent), .groups = "drop")

  annotColumns <- c(cf$table$fileName,
                    cf$table$sampleName,
                    cf$table$hkeysDepth(),
                    cf$table$fkeysDepth(),
                    cf$table$isotopeLabel)
  annotation <- data %>%
    dplyr::select(!!!syms(annotColumns) ) %>%
    distinct()

  res <- inner_join(annotation, pepStats, by = c(cf$table$sampleName, cf$table$hkeysDepth() ))
  res <-  if (nested) {res %>% group_by_at(cf$table$hkeysDepth()) %>% nest()} else {res}
  return(res)
}


# Summarize Intensities by Intensity or NAs ----
.rankProteinPrecursors <- function(data,
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
#' library(prolfqua)
#' bb <- prolfqua_data('data_spectronautDIA250_A')
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' res <- removeLarge_Q_Values(analysis, config)
#' #debug(rankPrecursorsByIntensity)
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
#' bb <- prolfqua_data('data_spectronautDIA250_A')
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' res <- removeLarge_Q_Values(analysis, config)
#' res <- rankPrecursorsByNAs(res,config)
#' colnames(res)
#' x <- res %>%
#'   dplyr::select(config$table$hierarchyKeys()[1], config$table$hierarchyKeys(TRUE)[1], "srm_NrNotNAs") %>%
#'   distinct() %>% dplyr::summarize(sum(srm_NrNotNAs)) %>% dplyr::pull()
#' stopifnot(sum(!is.na(res[[config$table$getWorkIntensity()[1]]])) == x)
#' res %>% dplyr::select(c(config$table$hierarchyKeys(),"srm_NrNotNAs"  ,"srm_NrNotNARank")) %>%
#'  distinct() %>%
#'  arrange(!!!syms(c(config$table$hierarchyKeys()[1],"srm_NrNotNARank")))
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
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua_data('data_spectronautDIA250_A')
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- analysis
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

  summaryPerPrecursorFiltered <- summaryPerPrecursor %>% dplyr::filter(.data$fraction > percent)
  summaryPerPrecursorFiltered <- summaryPerPrecursorFiltered %>%
    dplyr::select(c(table$hierarchyKeys())) %>% dplyr::distinct()
  stopifnot(all(colnames(summaryPerPrecursorFiltered) %in% table$hierarchyKeys()))
  res <- summaryPerPrecursorFiltered %>% left_join(pdata)
  return( dplyr::ungroup(res))
}


