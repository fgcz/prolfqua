# AnalysisParameters ----
#' Analysis parameters
#' @description
#' Analysis parameters
#' @keywords internal
#' @family configuration
#' @export
AnalysisParameters <- R6::R6Class(
  "AnalysisParameters",
  public = list(
    #' @field qVal_individual_threshold qValue threshold for sample
    qVal_individual_threshold  = 0.05,
    #' @field qVal_experiment_threshold qValue threshold for dataset
    qVal_experiment_threshold = 0.01,
    #' @field qVal_minNumber_below_experiment_threshold how many samples need to meet qVal_experiment_threshold
    qVal_minNumber_below_experiment_threshold = 3,
    #' @field min_nr_of_notNA minimum number of not NA's in all samples default 1
    min_nr_of_notNA = 1, # how many values per transition total
    #' @field min_nr_of_notNA_condition minimum number of not NA's in interaction
    min_nr_of_notNA_condition = 0, # how many not missing in condition
    #' @field min_peptides_protein minimum number of peptides per protein
    min_peptides_protein = 2
  )
)

