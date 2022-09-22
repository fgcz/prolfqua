#LFQData ----
#'
#' LFQData R6 class
#' @export
#' @family LFQData
#' @examples
#' istar <- old2new(prolfqua_data('data_ionstar')$filtered())
#'
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' #LFQData$debug("rename_response")
#' lfqdata <- LFQData$new(data, istar$config)
#' tmp <- lfqdata$to_wide()
#' stopifnot("data.frame" %in% class(tmp$data))
#' tmp <- lfqdata$to_wide(as.matrix = TRUE)
#' stopifnot("matrix" %in% class(tmp$data))
#' lfqdata$factors()
#' stopifnot(lfqdata$is_transformed()==FALSE)
#' lfqdata$summarize_hierarchy()
#' lfqdata$omit_NA()
#'
#' lfqdata$response()
#' lfqdata$rename_response("peptide.intensity")
#' lfqdata$response()
#' lfqdata$get_Plotter()$heatmap()
#' stopifnot("LFQData" %in% class(lfqdata$get_copy()))
#' stopifnot("LFQDataTransformer" %in% class(lfqdata$get_Transformer()))
#' stopifnot("LFQDataStats" %in% class(lfqdata$get_Stats()))
#' stopifnot("LFQDataSummariser" %in% class(lfqdata$get_Summariser()))
#' stopifnot("LFQDataPlotter" %in% class(lfqdata$get_Plotter()))
#' stopifnot("LFQDataWriter" %in% class(lfqdata$get_Writer()))
#' stopifnot("LFQDataAggregator" %in% class(lfqdata$get_Aggregator()))
#'
#'
LFQData <- R6::R6Class(
  "LFQData",

  public = list(
    #' @field config AnalysisConfiguration
    config = NULL,
    #' @field data data.frame or tibble matching AnalysisConfiguration.
    data = NULL,
    #' @field is_pep todo
    is_pep = FALSE,
    #' @field prefix e.g. "peptide_", "protein_", "compound_"
    prefix = "",
    #' @description
    #' initialize
    #' @param data data.frame
    #' @param config configuration
    #' @param is_pep todo
    #' @param prefix will be use as output prefix
    initialize = function(data, config, is_pep=TRUE, prefix = "ms_") {
      self$data <- data
      self$config <- config$clone(deep = TRUE)
      self$is_pep <- is_pep
      self$prefix <- prefix
    },
    #' @description
    #' get deep copy
    get_copy = function(){
      return(self$clone(deep = TRUE))
    },
    #' @description
    #' samples subset of data
    #' @param size size of subset default 100
    #' @param seed set seed
    get_sample = function(size = 100, seed = NULL){
      if(!is.null(seed)){ set.seed( seed )}
      subset <- prolfqua::sample_subset(size = size, self$data, self$config)
      return(LFQData$new(subset, self$config$clone(deep = TRUE)))
    },
    #' @description
    #' get subset of data
    #' @param x data frame with columns containing subject_Id
    get_subset = function(x){
      x <- select(x, all_of(self$subject_Id())) |> distinct()
      subset <- inner_join(x, self$data)
      return(LFQData$new(subset, self$config$clone(deep = TRUE)))
    },
    #' @description
    #' get subject ID columns
    subject_Id = function(){
      return(self$config$table$hkeysDepth())
    },
    #' @description
    #' is data trasfromed
    #' @param is_transformed logical
    #' @return logical
    is_transformed = function(is_transformed){
      if (missing(is_transformed)) {
        return(self$config$table$is_intensity_transformed)
      }else{
        self$config$table$is_intensity_transformed = is_transformed
      }
    },
    #' @description
    #' return name of intensity column
    #' @return name of intensity column
    intensity_column = function(){
      return(self$config$table$get_work_intensity())
    },
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    remove_small_intensities = function(threshold = 4){
      self$data <- prolfqua::remove_small_intensities( self$data, self$config, threshold = threshold )
      invisible(self)
    },
    #' @description
    #' remove proteins with less than X peptides
    #' @return self
    filter_proteins_by_peptide_count = function(){
      message("removing proteins with less than",
              self$config$parameter$min_peptides_protein,
              " peptpides")
      self$data <- prolfqua::filter_proteins_by_peptide_count(self$data, self$config)$data
      invisible(self)
    },
    #' @description
    #' Omit NA from intensities per hierarchy (e.g. protein or peptide), idea is to use it for normalization
    #' For instance if a peptide has a missing value in more then nrNA of the samples within a condition
    #' it will be removed
    #' @param nrNA number of NA values
    #' @param factorDepth you control for nrNA per condition or experiment etc. e.g. factorDepth = 0  then per experiment
    #' @return LFQData with NA omitted.
    omit_NA = function(nrNA = 0, factorDepth = FALSE){
      if (!factorDepth)
      {
        missing <- interaction_missing_stats(self$data, self$config)
      } else {
        cfg <- self$config$clone(deep = TRUE)
        cfg$table$factorDepth <- factorDepth
        missing <- interaction_missing_stats(self$data, cfg)
      }
      notNA <- missing$data |> dplyr::filter(nrNAs <= nrNA)
      sumN <- notNA |> group_by_at(self$config$table$hierarchyKeys()) |>
        summarise(n = n())
      notNA <- sumN |> dplyr::filter(n == max(n))

      notNA <- notNA |> dplyr::select(self$config$table$hierarchyKeys())
      notNAdata <- dplyr::inner_join( notNA, self$data) |> ungroup()
      return(LFQData$new(notNAdata, self$config$clone(deep = TRUE)))
    },
    #'
    #' @description
    #' some software is reporting NA's as 0, you must remove it from your data
    #' @param threshold default 4.
    #' @return self
    complete_cases = function(){
      self$data <- prolfqua::complete_cases(self$data, self$config)
      invisible(self)
    },
    #' @description
    #' converts the data to wide
    #' @param as.matrix return as data.frame or matrix
    #' @return list with data, annotation, and configuration
    to_wide = function(as.matrix = FALSE){
      wide <- prolfqua::tidy_to_wide_config(self$data, self$config, as.matrix = as.matrix)
      wide$config <- self$config$clone(deep = TRUE)
      return(wide)
    },
    #' @description
    #' Annotation table.
    #' @return data.frame
    factors = function(){
      prolfqua::table_factors(self$data, self$config)
    },
    #' @description
    #' name of response variable
    #' @return data.frame
    response = function(){
      self$config$table$get_work_intensity()
    },
    #' @description
    #' new name of response variable
    #' @param newname default Intensity
    rename_response = function(newname = "Intensity"){
      if((newname %in% colnames(self$data))){
        msg <- paste(newname, " already in data :", paste( colnames(self$data), collapse = " "), ".")
        logger::log_info(msg)
        logger::log_error("provide different name.")
      } else {
        old <- self$config$table$pop_work_intensity()
        self$config$table$set_work_intensity(newname)
        self$data <- self$data |> dplyr::rename(!!newname := !!sym(old))
      }
    },
    #' @description
    #' number of elements at each level
    hierarchy_counts = function(){
      prolfqua::hierarchy_counts(self$data, self$config)
    },
    #' @description
    #' e.g. number of peptides per protein etc
    #' @return data.frame
    summarize_hierarchy = function(){
      prolfqua::summarize_hierarchy(self$data, self$config)
    },
    #' @description
    #' get Plotter
    #' @return LFQDataPlotter
    get_Plotter = function(){
      return(LFQDataPlotter$new(self, self$prefix))
    },
    #' @description
    #' get Writer
    #' @param format array of formats to write to supported are xlsx, csv and html
    #' @return LFQDataPlotter
    get_Writer = function(format = "xlsx"){
      return(LFQDataWriter$new(self, self$prefix, format = format))
    },
    #' @description
    #' get Summariser
    #' @return LFQDataSummarizer
    get_Summariser = function(){
      return(LFQDataSummariser$new(self))
    },
    #' @description
    #' Get \code{\link{LFQDataStats}}. For more details see \code{\link{LFQDataStats}}.
    #' @param stats default interaction, computes statistics within interaction.
    #' @return LFQDataStats
    get_Stats = function(stats = c("everything","interaction", "all")){
      stats <- match.arg(stats)
      return(LFQDataStats$new(self, stats = stats))
    },
    #' @description
    #' get Stats
    #' @return LFQDataTransformer
    get_Transformer = function() {
      return(LFQDataTransformer$new(self))
    },
    #' @description
    #' get Aggregator
    #' @return LFQDataAggregator
    get_Aggregator = function() {
      return(LFQDataAggregator$new(self))
    },
    #' @description
    #' get difference of self with other if other is subset of self
    #'
    #' @details
    #' Use to compare filtering results obtained from self, e.g. which proteins and peptides were removed (other)
    #'
    #' @param other a filtered LFQData set
    #' @return LFQData
    #'
    filter_difference = function(other){
      diffdata <- prolfqua::filter_difference(self$data,other$data,self$config )
      res <- LFQData$new(diffdata , self$config$clone(deep = TRUE))
      return(res)
    }
  )
)





#' converts LFQData object to SummarizedExperiment
#'
#' For compatibility with Bioconductor
#' @param lfqdata LFQData object
#' @return SummarizedExperiment (bioconductor)
#' @family LFQData
#' @export
#'
LFQDataToSummarizedExperiment <- function(lfqdata){
  if (requireNamespace("SummarizedExperiment")) {
    wide <- lfqdata$to_wide(as.matrix = TRUE)
    ann <- data.frame(wide$annotation)
    rownames(ann) <- wide$annotation[[lfqdata$config$table$sampleName]]
    se <- SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(LFQ = wide$data), colData = ann)
    return(se)
  }
}
