
#' Add Annotation to a data.frame in long format
#' for an usage example see run_script lfq_mixed_model_inference
#' @family setup
#'
#' @param intensityData data imported using ``
#' @param inputAnnotation annotation
#' @param fileName column name to join on.
#' @export
#' @examples
#' protein_txt <- system_file("samples/maxquant_txt/tiny2.zip",package = "prolfqua")
#'
#' inputAnnotation <- system_file("samples/maxquant_txt/annotation_Ionstar2018_PXD003881.xlsx",package = "prolfqua")
#' startdata <- prolfqua::tidyMQ_ProteinGroups(protein_txt)
#' tmp <- add_annotation(startdata,inputAnnotation )
#' stopifnot(ncol(tmp) == ncol(startdata) + 3)
#'
add_annotation <- function(intensityData,
                           inputAnnotation,
                           fileName = "raw.file") {
  ## read the data
  {# add annotation
    if ( is.character(inputAnnotation) ) {
      annotation <- readxl::read_xlsx(inputAnnotation)
    } else {
      annotation <- inputAnnotation
    }
    noData <- annotation[!annotation[[fileName]] %in% intensityData[[fileName]],]
    if (nrow(noData)) {
      message("some files in annotation have no measurements")
      message(paste(noData, collapse = " "))
    }
    measSamples <- unique(intensityData[[fileName]])
    noAnnot <- measSamples[!measSamples %in% annotation[[fileName]] ]
    if (length(noAnnot) > 0 ) {
      message("some measured samples have no annotation!")
      message(paste(noAnnot,collapse = " "))
    }
    resPepProtAnnot <- inner_join(annotation, intensityData, by = fileName)
    ###  Setup analysis ####
  }
  return(resPepProtAnnot)
}


#' correlation preprocessing
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param minCorrelation correlation threshold default 0.7
#' @export
#' @keywords internal
#' @family workflows
#' @family deprecated
#' @examples
#'
#'
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' data$nr_peptide_Id_IN_protein_Id <- NULL
#'
#' config$parameter$min_nr_of_notNA  <- 3
#' #debug(workflow_correlation_preprocessing_protein_intensities)
#' runLong <- TRUE
#' if(runLong){
#' res <- workflow_correlation_preprocessing_protein_intensities(data,config)
#' names(res)
#' }
#'
#'
workflow_correlation_preprocessing_protein_intensities <- function(pdata, config, minCorrelation = 0.7){
  stat_input <- hierarchy_counts(pdata, config)

  data_NA_QVal <- filter_byQValue(pdata, config)
  stat_qval <- hierarchy_counts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- rank_by_NA(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal |> dplyr::filter(.data$srm_NrNotNAs > config$parameter$min_nr_of_notNA)
  if(nrow(data_NA_QVal) == 0){
    warning("no rows left after filtering for min_nr_of_notNA")
  }
  stat_min_nr_of_notNA <- hierarchy_counts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- filter_proteins_by_peptide_count(data_NA_QVal, config)$data

  stat_min_peptides_protein  <- hierarchy_counts(data_NA_QVal, config)

  # filter decorrelated.

  data_NA_QVal <- transform_work_intensity(data_NA_QVal, config, log2)
  data_NA_QVal <- mark_decorelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, .data$srm_decorelated == FALSE)

  stat_correlated  <- hierarchy_counts(keepCorrelated, config)

  # TODO check if you are not aggregating log transformed intensities
  # rank precursors by intensity
  keepCorrelated <- rank_peptide_by_intensity(keepCorrelated, config)
  qvalFiltImputed <- impute_correlationBased(keepCorrelated, config)
  mean_na <- function(x, name=FALSE){if(name){return("mean_na")};mean(x, na.rm = TRUE)}
  proteinIntensities <- aggregate_intensity_topN(qvalFiltImputed, config, .func = mean_na, N = 3)

  # collect stats
  stats <- list(stat_input = stat_input,
                stat_qval = stat_qval,
                stat_min_nr_of_notNA = stat_min_nr_of_notNA,
                stat_min_peptides_protein = stat_min_peptides_protein,
                stat_correlated = stat_correlated
  )
  x <- dplyr::bind_rows(stats)
  stats <- add_column(x, processing = names(stats),.before = 1)


  return(list(data = proteinIntensities$data, stats = stats, newconfig = proteinIntensities$newconfig))
}

#' Apply correlation filtering and impute missing values
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param minCorrelation correlation threshold default 0.7
#' @keywords internal
#' @family workflows
#' @export
#' @examples
#'
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' data$nr_peptide_Id_IN_protein_Id <- NULL
#'
#' config$parameter$min_nr_of_notNA  <- 3
#' res <- workflow_corr_filter_impute(data,config)
#'
workflow_corr_filter_impute <- function(pdata, config, minCorrelation =0.6){
  stat_input <- hierarchy_counts(pdata, config)

  data_NA_QVal <- filter_byQValue(pdata, config)
  stat_qval <- hierarchy_counts(data_NA_QVal, config)

  # remove transitions with large numbers of NA's
  data_NA_QVal <- rank_by_NA(data_NA_QVal, config)
  data_NA_QVal <- data_NA_QVal |> dplyr::filter(.data$srm_NrNotNAs > config$parameter$min_nr_of_notNA)
  stat_min_nr_of_notNA <- hierarchy_counts(data_NA_QVal, config)

  # remove single hit wonders
  data_NA_QVal <- filter_proteins_by_peptide_count(data_NA_QVal, config)$data

  stat_min_peptides_protein  <- hierarchy_counts(data_NA_QVal, config)

  # filter decorrelated.
  data_NA_QVal <- transform_work_intensity(data_NA_QVal, config, log2)
  data_NA_QVal <- mark_decorelated(data_NA_QVal, config, minCorrelation = minCorrelation)
  keepCorrelated <- dplyr::filter(data_NA_QVal, .data$srm_decorelated == FALSE)

  stat_correlated  <- hierarchy_counts(keepCorrelated, config)
  keepCorrelated <- rank_peptide_by_intensity(keepCorrelated, config)
  qvalFiltImputed <- impute_correlationBased(keepCorrelated, config)


  stats <- list(stat_input = stat_input,
                stat_qval = stat_qval,
                stat_min_nr_of_notNA = stat_min_nr_of_notNA,
                stat_min_peptides_protein = stat_min_peptides_protein,
                stat_correlated = stat_correlated
  )
  x <- dplyr::bind_rows(stats)
  stats <- add_column(x, processing = names(stats),.before = 1)
  return(qvalFiltImputed)
}

#' filter QVAlues and NA's and factor information
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param minCorrelation correlation threshold default 0.7
#' @export
#' @keywords internal
#' @family workflows
#' @family deprecated
#' @examples
#'
#' rm(list=ls())
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#' data$nr_peptide_Id_IN_protein_Id <- NULL
#'
#' hierarchy_counts(data, config)
#' tmp <- workflow_DIA_NA_preprocessing(data, config)
#' hierarchy_counts(tmp$data, config)
#' tmp <- workflow_DIA_NA_preprocessing(data, config, percent=70)
#' hierarchy_counts(tmp$data, config)
#' stopifnot(FALSE == (dplyr::is_grouped_df(tmp$data)))
#'
workflow_DIA_NA_preprocessing <- function(pdata,
                                          config,
                                          percent = 60,
                                          hierarchy_level = 2,
                                          factor_level = 1,
                                          min_peptides_protein = config$parameter$min_peptides_protein)
{
  stat_input <- hierarchy_counts(pdata, config)

  data_NA_QVal <- filter_byQValue(pdata, config)
  stat_qval <- hierarchy_counts(data_NA_QVal, config)

  resNACondition <- filter_factor_levels_by_missing(data_NA_QVal,
                                                    config,
                                                    percent = percent)

  stat_naFilter <- hierarchy_counts(resNACondition, config)
  protID <- summarize_hierarchy(resNACondition,config) |>
    dplyr::filter(!!sym(paste0(config$table$hierarchy_keys()[hierarchy_level],"_n"))
                  >= min_peptides_protein)

  data_NA_QVal_condition <- protID |>
    dplyr::select(config$table$hierarchy_keys()[1]) |>
    dplyr::inner_join(resNACondition)

  # Complete cases
  data_NA_QVal_condition <- complete_cases( data_NA_QVal_condition , config)
  stat_peptidFitler <- hierarchy_counts(data_NA_QVal_condition, config)
  stats = list(stat_input = stat_input,
               stat_qval = stat_qval,
               stat_naFilter = stat_naFilter,
               stat_peptidFitler = stat_peptidFitler
  )
  return(list(data = data_NA_QVal_condition, stats = stats))
}

