# Direct intensity manipulation ----

#' Remove rows when qVal_individual_threshold is exceeded
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
#' res <- remove_large_QValues(analysis, config)
remove_large_QValues <- function(pdata, config, qValThreshold = config$parameter$qVal_individual_threshold){
  pdata <- pdata |>
    dplyr::filter(!!sym(config$table$ident_qValue) < qValThreshold)
  return(pdata)
}

#' Remove rows when intensity lower then threshold
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
#' config$table$get_response()
#'
#' res1 <- remove_small_intensities(analysis, config, threshold=1 )
#' res1000 <- remove_small_intensities(analysis, config, threshold=1000 )
#' stopifnot(nrow(res1) >  nrow(res1000))
#'
remove_small_intensities <- function(pdata, config, threshold = 1){
  resData <- pdata |> dplyr::filter(!!sym(config$table$get_response()) >= threshold)
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
#'
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#' x <- transform_work_intensity(analysis, config, .func = log2)
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
    newcol <- paste(.funcname, config$table$get_response(), sep = "_")
  }else{
    newcol <- intesityNewName
  }

  #pdata <- pdata |> dplyr::mutate_at(config$table$get_response(),
  #                                    .funs = funs(!!sym(newcol) := .func(.)))
  pdata <- pdata |> dplyr::mutate(!!sym(newcol) := .func(!!sym(config$table$get_response())))

  config$table$set_response(newcol)
  message("Column added : ", newcol)
  config$table$is_response_transformed = TRUE

  if (deep) {
    return( list(data = pdata, config = config) )
  } else {
    return( pdata)
  }
}



# Summarize Q Values ----
#'
#' Compute QValue summaries for each precursor or peptide or protein
#'
#'
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
#' res <- summarise_QValues(analysis, config)
#' stopifnot(c("srm_QValueMin", "srm_QValueNR") %in% colnames(res))
#' hist(unique(res$srm_QValueMin))
#' hist(unique(res$srm_QValueNR))
#'
summarise_QValues <- function(pdata,
                              config
){
  QValueMin <- "srm_QValueMin"
  QValueNR <- "srm_QValueNR"

  precursorIDs <- config$table$hierarchy_keys()
  fileName <- config$table$fileName
  QValue  <- config$table$ident_qValue
  qVal_minNumber_below_experiment_threshold <- config$parameter$qVal_minNumber_below_experiment_threshold
  qVal_experiment_threshold <- config$parameter$qVal_experiment_threshold

  nthbestQValue <-  function(x,qVal_minNumber_below_experiment_threshold){sort(x)[qVal_minNumber_below_experiment_threshold]}
  npass <-  function(x,thresh = qVal_experiment_threshold){sum(x < thresh)}

  qValueSummaries <- pdata |>
    dplyr::select(!!!syms(c(fileName, precursorIDs, config$table$ident_qValue))) |>
    dplyr::group_by_at(precursorIDs) |>
    dplyr::summarise_at(  c( QValue ), .funs = funs(!!QValueMin := nthbestQValue(.,qVal_minNumber_below_experiment_threshold ),
                                                    !!QValueNR  := npass(., qVal_experiment_threshold)
    ))
  pdata <- dplyr::inner_join(pdata, qValueSummaries, by = c(precursorIDs))
  message("Columns added :",QValueMin, QValueNR)
  return(pdata)
}

#' Filter data by individual and experiment wide thresholds
#'
#' employs parameters ident_qValue, qVal_minNumber_below_experiment_threshold,
#' qVal_individual_threshold and qVal_experiment_threshold
#' @export
#' @keywords internal
#' @examples
#'
#' dd <- prolfqua_data('data_spectronautDIA250_A')
#' config <- dd$config_f()
#' analysis <- dd$analysis(dd$data,config)
#'
#' summarize_hierarchy(analysis, config)
#' res <- filter_byQValue(analysis, config)
#' summarize_hierarchy(res, config)
#'
filter_byQValue <- function(pdata, config){
  data_NA <- remove_large_QValues(pdata, config)
  data_NA <- summarise_QValues(data_NA, config)
  data_NA_QVal <- data_NA |>
    dplyr::filter( !!sym("srm_QValueMin") < config$parameter$qVal_experiment_threshold )
}


# Intensities to wide ----

.ExtractMatrix <- function(x, sep = "~lfq~"){
  idx <- sapply(x,is.numeric)
  xmat <- as.matrix(x[,idx])
  idcols <- x |> dplyr::select(which(!idx == TRUE))
  if (ncol(idcols) > 0) {
    rownames(xmat) <- x |> dplyr::select(which(!idx == TRUE)) |>
      tidyr::unite(x, sep = sep) |> dplyr::pull(x)
  }
  xmat
}



#' Transform tidy table into a table with a column of responses for each sample
#'
#' @export
#' @keywords internal
tidy_to_wide <- function(data,
                         rowIDs ,
                         columnLabels ,
                         value
){
  wide <- data |>
    dplyr::select_at(c(rowIDs, columnLabels, value  ))

  wide_spread <- wide |>
    tidyr::spread(key = columnLabels, value = value)
  # wide_pivot <- wide |>
  #  tidyr::pivot_wider( names_from = all_of(columnLabels) , values_from =  all_of(value) )

  #message("I AM BEING USED: spread")
  return(wide_spread)
}

#' transform long to wide
#' @export
#' @keywords internal
#' @return list with data, rowdata, and annotation (colData)
#' @examples
#'
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' config <- dd$config
#' data <- dd$data
#' res <- tidy_to_wide_config(data, config)
#' testthat::expect_equal(nrow(res$rowdata), nrow(res$data))
#' testthat::expect_equal(ncol(res$data) - ncol(res$rowdata) , nrow(res$annotation))
#' res <- tidy_to_wide_config(data, config, as.matrix = TRUE)
#' stopifnot(all(dim(res$data) == c(28,  12)))
#' stopifnot(all(dim(res$annotation) == c(12,  4)))
#' stopifnot(all(dim(res$rowdata) == c(28, 3)))
#'
#' res <- scale(res$data)
#' tidy_to_wide_config(data, config,  value = config$table$nr_children)
#'
#'
#' xt <- prolfqua::LFQData$new(dd$data, dd$config)
#' xt$data$nr_children
#' #xt$config$table$is_response_transformed <- TRUE
#' res <- xt$get_Aggregator()
#' x <- res$medpolish()
#' towide <- tidy_to_wide_config(x$data, x$config,  value = x$config$table$nr_children)
#'
#' dd <- prolfqua::sim_lfq_data_protein_config()
#' dd$config$table$nr_children
#' dd$data
#' xt <- tidy_to_wide_config(dd$data, dd$config,  value = dd$config$table$nr_children)
#' xt$data
#'
tidy_to_wide_config <- function(data, config,
                                as.matrix = FALSE,
                                fileName = FALSE,
                                sep="~lfq~",
                                value = config$table$get_response()
){
  if (fileName) {
    newcolname <- config$table$fileName
  }else{
    newcolname <- config$table$sampleName
  }

  ids <- dplyr::select(data,
                       all_of(c( config$table$sampleName, config$table$fileName, config$table$factor_keys(), config$table$isotopeLabel))) |>
    dplyr::distinct() |> dplyr::arrange_at(newcolname)

  res <- tidy_to_wide( data, c(config$table$hierarchy_keys(),config$table$isotopeLabel) ,
                       newcolname,
                       value = value )
  rowdata <- res |> dplyr::select(all_of(c(config$table$hierarchy_keys(),config$table$isotopeLabel)))
  if (as.matrix) {
    resMat <- as.matrix(dplyr::select(res,-dplyr::one_of(c(config$table$hierarchy_keys(),config$table$isotopeLabel))))
    names <- rowdata |>
      tidyr::unite("newID", !!!dplyr::syms(c(config$table$hierarchy_keys(), config$table$isotopeLabel)), sep = sep) |> dplyr::pull("newID")
    rownames(resMat) <- names
    res <- resMat
  }
  return(list(data = res, annotation = ids, rowdata = rowdata))
}

#' Takes matrix of responses and converts into tibble
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
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' data <- dd$data
#' conf <- dd$config
#' res <- tidy_to_wide_config(data, conf, as.matrix = TRUE)
#'
#' res <- scale(res$data)
#' xx <- response_matrix_as_tibble(res,"srm_intensityScaled", conf)
#' xx <- response_matrix_as_tibble(res,"srm_intensityScaled", conf, data)
#' conf$table$get_response() == "srm_intensityScaled"
#'
response_matrix_as_tibble <- function(pdata, value, config, data = NULL, sep = "~lfq~"){
  pdata <- dplyr::bind_cols(
    tibble::tibble("row.names" := rownames(pdata)),
    tibble::as_tibble(pdata)
  )
  pdata <- tidyr::gather(pdata,key = !!config$table$sampleName, value = !!value, 2:ncol(pdata))
  pdata <- tidyr::separate(pdata, "row.names",  config$table$hierarchy_keys(), sep = sep)
  if (!is.null(data)) {
    pdata <- dplyr::inner_join(data, pdata)
    config$table$set_response(value)
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
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' conf <- bb$config
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
  data <- tidy_to_wide_config(data, config, as.matrix = TRUE)$data
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


#' Apply function requiring a matrix to tidy table
#'
#' @param data data.frame
#' @param config AnalysisConfiguration
#' @param .func function
#' @param .funcname name of function (used for creating new column)
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' data <- bb$data
#' conf <- bb$config
#' res <- apply_to_response_matrix(data, conf, .func = base::scale)
#'
#' stopifnot("abundance_base..scale" %in% colnames(res))
#' stopifnot("abundance_base..scale" == conf$table$get_response())
#' conf <- bb$config$clone(deep=TRUE)
#' conf$table$workIntensity <- "abundance"
#' res <- apply_to_response_matrix(data, conf$clone(deep=TRUE), .func = robust_scale)
#'
#' # Normalize data using the vsn method from bioconductor
#'
#' if( require("vsn")){
#'  res <- apply_to_response_matrix(data, conf$clone(deep=TRUE), .func = vsn::justvsn)
#' }
#'
apply_to_response_matrix <- function(data, config, .func, .funcname = NULL){
  .call <- as.list( match.call() )
  .funcname <- if (is.null(.funcname)) { deparse(.call$.func) } else {.funcname}
  colname <- make.names( paste( config$table$get_response(), .funcname, sep = "_"))
  mat <- tidy_to_wide_config(data, config, as.matrix = TRUE)$data
  mat <- .func(mat)
  data <- response_matrix_as_tibble(mat, colname, config, data)
  return(data)
}

#' Scale data using a subset of the data
#'
#' this should reduce the overall variance.
#'
#' @export
#' @keywords internal
#' @param data the whole dataset
#' @param subset a subset of the dataset
#' @param config configuration
#' @param perserveMean default FASE - sets mean to zero
#' @param get_scales return a list of transformed data and the scaling parameters
#' @family preprocessing
#' @examples
#'
#'
#'
#' bb <-sim_lfq_data_peptide_config(Nprot = 100)
#' conf <- bb$config$clone(deep=TRUE)
#' sample_analysis <- bb$data
#'
#' res <- transform_work_intensity(sample_analysis, conf, log2)
#' s1 <- get_robscales(res, conf)
#' res <- scale_with_subset(res, res, conf)
#' s2 <- get_robscales(res$data, conf)
#' stopifnot(abs(mean(s1$mads) - mean(s2$mads)) < 1e-6)
scale_with_subset <- function(data, subset, config, preserveMean = FALSE, get_scales = TRUE){

  colname <- make.names( paste( config$table$get_response(), "subset_scaled", sep = "_"))
  subset <- tidy_to_wide_config(subset, config, as.matrix = TRUE)$data


  scales <- .get_robscales(subset)
  mat <- tidy_to_wide_config(data, config, as.matrix = TRUE)$data
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
  data <- response_matrix_as_tibble(mat, colname, config, data)
  if (get_scales) {
    return(list(data = data, scales = scales))
  } else {
    return(data)
  }
}

#' Scale using a subset of the data, within factor levels (e.g. use for pulldown data)
#'
#' This method reduces the variance within the group.
#'
#' @export
#' @keywords internal
#' @param data tibble with data
#' @param subset tibble with subset of the data which will be used to derive scaling parameters
#' @param config configuration
#' @param preserveMean default FALSE then set mean to 0
#' @family preprocessing
#' @examples
#'
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' conf <- bb$config$clone(deep=TRUE)
#' sample_analysis <- bb$data
#'
#' res <- transform_work_intensity(sample_analysis, conf, log2)
#' res <- scale_with_subset_by_factors(res, res, conf)
#'
#'
scale_with_subset_by_factors <-  function(data, subset, config, preserveMean = FALSE){
  config <- config$clone(deep = TRUE)
  dl <- group_by(data, !!!syms(config$table$factor_keys_depth())) |> nest()
  sl <- group_by(subset, !!!syms(config$table$factor_keys_depth())) |> nest()
  cf <- config$clone(deep = TRUE)
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
  config$table$set_response(cf$table$get_response())
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
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' xx <- normalize_log2_robscale(bb$data, bb$config)
#' xx$config$table$workIntensity
#'
normalize_log2_robscale <- function(pdata, config){
  pepConfig <- config$clone(deep = TRUE)
  pepIntensityNormalized <- transform_work_intensity(pdata, pepConfig, log2)
  pepConfig$table$is_response_transformed = TRUE

  pepIntensityNormalized <- apply_to_response_matrix(pepIntensityNormalized,
                                                     pepConfig,
                                                     .func = robust_scale)

  pepIntensityNormalized <- pepIntensityNormalized |>
    dplyr::rename(transformedIntensity = pepConfig$table$get_response())
  pepConfig$table$pop_response()
  pepConfig$table$set_response("transformedIntensity")

  return(list(data = pepIntensityNormalized, config = pepConfig))
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
  if (!c_name %in% colnames(data) ) {
    data$c_name <- NULL
  }
  tmp <- data |>
    dplyr::select_at(c(levelA, levelB)) |>
    dplyr::distinct() |>
    dplyr::group_by_at(levelA) |>
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
#'
#' bb <- sim_lfq_data_peptide_config(Nprot = 100)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' hierarchy <- config$table$hierarchy_keys()
#' res <- nr_B_in_A(data, config)
#'
#' res$data |>
#'   dplyr::select_at(c(config$table$hierarchy_keys_depth(),  res$name)) |>
#'   dplyr::distinct() |>
#'   dplyr::pull() |> table()
#'
#'
#' bb <- prolfqua::prolfqua_data('data_skylineSRM_HL_A')
#' config <- old2new(bb$config_f())
#' data <- bb$data
#' data$Area[data$Area == 0] <- NA
#' analysis <- setup_analysis(data, config)
#'
#' resDataStart <- bb$analysis(bb$data, config)
#'
#' nr_B_in_A(resDataStart, config)
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#' config$table$hierarchyDepth <- 2
#' nr_B_in_A(resDataStart, config, merge = FALSE)
#'
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
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' nr_B_in_A_per_sample(bb$data, bb$config, nested =FALSE)
#' bb <- prolfqua::sim_lfq_data_protein_config()
#' nr_B_in_A_per_sample(bb$data, bb$config, nested=FALSE)
#'
nr_B_in_A_per_sample <- function(data, config, nested = TRUE){
  #TODO wew check for deprecation since not used.
  cf <- config

  levelA <- cf$table$hierarchy_keys_depth()
  levelB <- cf$table$hierarchy_keys()[length(levelA) + 1]
  if (is.na(levelB)) {
    warning("here is no B in A")
  }
  data <- prolfqua::complete_cases(data, cf)
  data <- data |>
    dplyr::mutate(presentabsent = case_when(!is.na(!!sym(cf$table$get_response())) ~ 1,
                                            TRUE ~ 0))
  pepStats <- data |> group_by_at(c(cf$table$hierarchy_keys_depth(), cf$table$sampleName)) |>
    summarize(nrPep = n(), present = sum(.data$presentabsent), .groups = "drop")

  annotColumns <- c(cf$table$fileName,
                    cf$table$sampleName,
                    cf$table$hierarchy_keys_depth(),
                    cf$table$factor_keys_depth(),
                    cf$table$isotopeLabel)
  annotation <- data |>
    dplyr::select(!!!syms(annotColumns) ) |>
    distinct()

  res <- inner_join(annotation, pepStats, by = c(cf$table$sampleName, cf$table$hierarchy_keys_depth() ))
  res <-  if (nested) {res |> group_by_at(cf$table$hierarchy_keys_depth()) |> nest()} else {res}
  return(res)
}



#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @export
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' dd$data <- na.omit(dd$data)
#' xd <- nr_obs_sample(dd$data, dd$config)
#' xd$nr_children |> table()
#' # xd |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_children)
#'
#' dp <- prolfqua::sim_lfq_data_protein_config()
#' xp <- nr_obs_sample(dp$data, dp$config)
#' # xp
#' # xp |> pivot_wider(id_cols = protein_Id, names_from = sample, values_from = nr_peptides)
#' xp$nr_peptides |> table()
#'
nr_obs_sample <- function(data, config, new_child = config$table$nr_children){
  data <- na.omit(data)
  nr_children <- data |>
    group_by(!!!rlang::syms(c(config$table$hierarchy_keys_depth(), config$table$fileName))) |>
    summarize(!!new_child := sum(!!sym(config$table$nr_children), na.rm = TRUE), .groups = "drop")
  return(nr_children)
}

#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @export
#' @param data tidy data
#' @param config prolfqua config
#' @param from_children TRUE compute from nr_children column, FALSE - do not use nr_children column
#' @param name_nr_child how to name column
#' @examples
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#'
#'
#' xd <- nr_obs_experiment(dd$data, dd$config, from_children = FALSE)
#' xd <- nr_obs_experiment(dd$data, dd$config, from_children = TRUE)
#' stopifnot(min(xd$nr_child_exp) == 1)
#' dp <- prolfqua::sim_lfq_data_protein_config()
#' nr_obs_experiment(dp$data, dp$config)
#' nr_obs_experiment(dp$data, dp$config, from_children = FALSE)
#'
#'
#' dd <- prolfqua::sim_lfq_data_peptide_config()
#' dd$config$table$hierarchyDepth <- 2
#' xpep <- nr_obs_experiment(dd$data,dd$config)
#' stopifnot(all(xpep$nr_child_exp == 1))
#' xpep <- nr_obs_experiment(dd$data,dd$config, from_children = FALSE)
#' stopifnot(all(xpep$nr_child_exp == 1))
#'
nr_obs_experiment <- function(data, config,
                              from_children = TRUE,
                              name_nr_child = "nr_child_exp"){
  tb <- config$table
  hkeys <- tb$hierarchy_keys()
  hkeysd <- tb$hierarchy_keys_depth()
  nr_children <- tb$nr_children
    if (from_children ) {

    xz <- nr_obs_sample(data,config)
    xz <- xz |> group_by(!!!syms(hkeysd)) |>
      dplyr::summarize(!!name_nr_child := max(!!sym(nr_children)), .groups = "drop")
    return(xz)
  } else {
    xq <- data |> dplyr::select(hkeys) |>
      distinct() |>
      dplyr::group_by(!!!syms(hkeysd)) |>
      dplyr::summarize(!!name_nr_child := dplyr::n(), .groups = "drop")
    return(xq)
  }
}


# Summarize Intensities by Intensity or NAs ----
.rankProteinPrecursors <- function(data,
                                   config,
                                   column = config$table$get_response(),
                                   fun = function(x){ sum(x, na.rm = TRUE)},
                                   summaryColumn = "srm_meanInt",
                                   rankColumn = "srm_meanIntRank",
                                   rankFunction = function(x){ min_rank(desc(x)) }
){
  table <- config$table

  summaryPerPrecursor <- data |>
    dplyr::group_by(!!!syms(table$hierarchy_keys())) |>
    dplyr::summarize(!!summaryColumn := fun(!!sym(column)))

  groupedByProtein <- summaryPerPrecursor |>
    dplyr::arrange(!!sym( table$hierarchy_keys()[1])) |>
    dplyr::group_by(!!sym( table$hierarchy_keys()[1]))
  rankedBySummary <- groupedByProtein |>
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
#'
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' res <- remove_large_QValues(bb$data, bb$config)
#' res <- rank_peptide_by_intensity(res,bb$config)
#' X <-res |> dplyr::select(c(bb$config$table$hierarchy_keys(),
#'  srm_meanInt, srm_meanIntRank)) |> dplyr::distinct()
#' X |> dplyr::arrange(!!!rlang::syms(c(bb$config$table$hierarchy_keys()[1], "srm_meanIntRank"  )))
rank_peptide_by_intensity <- function(pdata, config){
  summaryColumn <- "srm_meanInt"
  rankColumn <- "srm_meanIntRank"
  pdata <- .rankProteinPrecursors(pdata, config, column = config$table$get_response(),
                                  fun = function(x){ mean(x, na.rm = TRUE)},
                                  summaryColumn = summaryColumn,
                                  rankColumn = rankColumn,
                                  rankFunction = function(x){min_rank(desc(x))}
  )

  message("Columns added : ", summaryColumn, " ",  rankColumn)
  return(pdata)
}

# Summarise NAs on lowest hierarchy ----

#' Ranks peptides/precursors of a protein by NAs (adds new column .NARank)
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @return data.frame
#' @export
#' @keywords internal
#' @examples
#'
#' bb <- prolfqua::prolfqua_data('data_spectronautDIA250_A')
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' res <- remove_large_QValues(analysis, config)
#' res <- rank_by_NA(res,config)
#' colnames(res)
#' x <- res |>
#'   dplyr::select(config$table$hierarchy_keys()[1], config$table$hierarchy_keys(TRUE)[1], "srm_NrNotNAs") |>
#'   dplyr::distinct() |> dplyr::summarize(sum(srm_NrNotNAs)) |> dplyr::pull()
#' stopifnot(sum(!is.na(res[[config$table$get_response()[1]]])) == x)
#' res |> dplyr::select(c(config$table$hierarchy_keys(),"srm_NrNotNAs"  ,"srm_NrNotNARank")) |>
#'  dplyr::distinct() |>
#'  dplyr::arrange(!!!rlang::syms(c(config$table$hierarchy_keys()[1],"srm_NrNotNARank")))
rank_by_NA <- function(pdata, config){
  summaryColumn <- "srm_NrNotNAs"
  rankColumn <- "srm_NrNotNARank"
  pdata <- .rankProteinPrecursors(pdata, config,
                                  column = config$table$get_response(),
                                  fun = function(x){sum(!is.na(x))},
                                  summaryColumn = summaryColumn,
                                  rankColumn = rankColumn,
                                  rankFunction = function(x){min_rank(desc(x))}
  )
  message("Columns added : ", summaryColumn, " ",  rankColumn)
  return(pdata)
}


