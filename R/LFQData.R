#LFQData ----
#'
#' LFQData R6 class
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' #LFQData$debug("omit_NA")
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' lfqdata$filter_proteins_by_peptide_count()
#' tmp <- lfqdata$to_wide()
#' testthat::expect_equal(nrow(tmp$data) , nrow(tmp$rowdata))
#' testthat::expect_equal(ncol(tmp$data) , nrow(tmp$annotation) + ncol(tmp$rowdata))
#'
#' stopifnot("data.frame" %in% class(tmp$data))
#' tmp <- lfqdata$to_wide(as.matrix = TRUE)
#' stopifnot("matrix" %in% class(tmp$data))
#' stopifnot(lfqdata$is_transformed()==FALSE)
#' lfqdata$summarize_hierarchy()
#'
#' # filter for missing values
#'
#' f1 <- lfqdata$omit_NA(nrNA = 0)
#' stopifnot(f1$hierarchy_counts() <= lfqdata$hierarchy_counts())
#'
#' f2 <- lfqdata$omit_NA(factorDepth = 0)
#' stopifnot(f2$hierarchy_counts() <= lfqdata$hierarchy_counts())
#'
#' lfqdata$response()
#' lfqdata$rename_response("peptide.intensity")
#' lfqdata$response()
#' stopifnot("LFQData" %in% class(lfqdata$get_copy()))
#' stopifnot("LFQDataTransformer" %in% class(lfqdata$get_Transformer()))
#' stopifnot("LFQDataStats" %in% class(lfqdata$get_Stats()))
#' stopifnot("LFQDataSummariser" %in% class(lfqdata$get_Summariser()))
#' stopifnot("LFQDataPlotter" %in% class(lfqdata$get_Plotter()))
#' stopifnot("LFQDataWriter" %in% class(lfqdata$get_Writer()))
#' stopifnot("LFQDataAggregator" %in% class(lfqdata$get_Aggregator()))
#'
#' lfqdata2 <- lfqdata$get_copy()
#' lfqdata2$data <- lfqdata2$data[1:100,]
#' res <- lfqdata$filter_difference(lfqdata2)
#' stopifnot(nrow(res$data) == nrow(lfqdata$data) - 100)
#'
#' tmp <- lfqdata$get_sample(5, seed = 4)
#' stopifnot(nrow(tmp$hierarchy()) == 5)
#'
#' lw <- lfqdata$get_Writer()
#' stopifnot(names(lw$get_wide()) %in% c("data", "annotation"))
#'
#' stopifnot("data.frame" %in% class(lw$get_long()))
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
    #' @param setup is data setup needed, default = FALSE, if TRUE, calls \code{\link{setup_analysis}} on data first.
    initialize = function(data, config, is_pep=TRUE, prefix = "ms_", setup = FALSE) {
      self$data <- if (setup) {setup_analysis(data, config)} else {data}
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
      if (!is.null(seed)) {  set.seed( seed ) }
      subset <- prolfqua::sample_subset(size = size, self$data, self$config)
      return(LFQData$new(subset, self$config$clone(deep = TRUE)))
    },
    #' @description
    #' get subset of data
    #' @param x data frame with columns containing subject_Id
    get_subset = function(x){
      x <- select(x, any_of(self$subject_Id())) |> distinct()
      subset <- inner_join(x, self$data)
      return(LFQData$new(subset, self$config$clone(deep = TRUE)))
    },
    #' @description
    #' get subject ID columns
    subject_Id = function(){
      return(self$config$table$hierarchy_keys_depth())
    },
    #' @description
    #' is data trasfromed
    #' @param is_transformed logical
    #' @return logical
    is_transformed = function(is_transformed){
      if (missing(is_transformed)) {
        return(self$config$table$is_response_transformed)
      }else{
        self$config$table$is_response_transformed = is_transformed
      }
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
      message("removing proteins with less than: ",
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
    #'
    omit_NA = function(nrNA = 0, factorDepth = NULL){
      if (is.null(factorDepth)) {
        missing <- prolfqua::summarize_stats_factors(self$data, self$config)
      } else {
        if (factorDepth >= 1) {
          cfg <- self$config$clone(deep = TRUE)
          cfg$table$factorDepth <- factorDepth
          missing <- prolfqua::summarize_stats_factors(self$data, cfg)
        } else{
          missing <- prolfqua::summarize_stats_all(self$data, self$config)
        }
      }
      notNA <- missing |> dplyr::filter(nrNAs <= nrNA)
      sumN <- notNA |> group_by_at(self$config$table$hierarchy_keys()) |>
        summarise(n = n())
      notNA <- sumN |> dplyr::filter(n == max(n))

      notNA <- notNA |> dplyr::select(self$config$table$hierarchy_keys())
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
    #' @param value either response or nr chidren
    #' @return list with data, annotation, and configuration
    to_wide = function(as.matrix = FALSE, value = c("response", "nr_children")){
      value <- match.arg(value)
      if (value == "response") {
        wide <- prolfqua::tidy_to_wide_config(self$data, self$config, as.matrix = as.matrix)
      } else {
        wide <- prolfqua::tidy_to_wide_config(
          self$data, self$config,
          as.matrix = as.matrix,
          value = self$config$table$nr_children)
      }
      wide$config <- self$config$clone(deep = TRUE)
      return(wide)
    },
    #' @description
    #' Annotation table
    #' @return data.frame
    factors = function(){
      prolfqua::table_factors(self$data, self$config)
    },
    #' @description
    #' Hierarchy table
    hierarchy = function(){
      hk <- self$config$table$hierarchy_keys_depth()
      hkdf <- self$data |> select(all_of(hk)) |> distinct()
      return(hkdf)
    },
    #' @description
    #' name of response variable
    #' @return data.frame
    response = function(){
      self$config$table$get_response()
    },
    #' @description
    #' new name of response variable
    #' @param newname default Intensity
    rename_response = function(newname = "Intensity"){
      if((newname %in% colnames(self$data))){
        msg <- paste(newname, " already in data :", paste( colnames(self$data), collapse = " "), ".")
        message(msg)
      } else {
        old <- self$config$table$pop_response()
        self$config$table$set_response(newname)
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
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' istar$config <- (istar$config)
#' data <- istar$data
#' lfqdata <- LFQData$new(data, istar$config)
#' lfqdata$to_wide()
#' if(require("SummarizedExperiment")){
#'    tmp <- LFQDataToSummarizedExperiment(lfqdata)
#' }
#'
LFQDataToSummarizedExperiment <- function(lfqdata){
  if (requireNamespace("SummarizedExperiment")) {
    wide <- lfqdata$to_wide(as.matrix = TRUE)
    nr_children <- lfqdata$to_wide(as.matrix = TRUE, value = "nr_children")
    ann <- data.frame(wide$annotation)
    rownames(ann) <- wide$annotation[[lfqdata$config$table$sampleName]]
    se <- SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(
      LFQ = wide$data,
      nr_children = nr_children$data),
      colData = ann,
      rowData = wide$rowdata)
    return(se)
  }
}
