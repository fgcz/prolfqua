# Helper functions -----

#' Keep only those proteins with 2 IDENTIFIED peptides
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return list with data.frame (data) and name of new column (name)
#' @export
#' @keywords internal
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#'
#' istar <- LFQServiceData::ionstar$Pep()
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' filterPep <- LFQService:::filter_proteins_by_peptide_count( istar$data ,  istar$config )
#'  x <- LFQService::summarize_hierarchy(filterPep$data , istar$config)
#' stopifnot(x$peptide_Id_n >= istar$config$parameter$min_peptides_protein)
#'
#'
filter_proteins_by_peptide_count <-
  function(pdata,
           config){

    # remove single hit wonders
    pdata <- nr_B_in_A(pdata,config)
    res <- dplyr::filter(pdata$data, !!sym(pdata$name) >= config$parameter$min_peptides_protein )

    return(list(data = res, name = pdata$name))
  }

#' normalize data by log2 and robust scaling
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return list with data.frame (data) and updated config (config)
#' @export
#' @keywords internal
#' @examples
#'
#' istar <- LFQServiceData::ionstar$Pep()
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' filterPep <- LFQService:::filter_proteins_by_peptide_count( istar$data ,  istar$config )
#' xx <- normalize_log2_robscale(istar$data, istar$config)
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


#' get the difference of two dataset where one is a subset of the other.
#' @param x data.frame
#' @param y data.frame
#' @param config AnlysisConfiguration
#' @return data.frame
#' @export
#' @keywords internal
#'
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#'
#' istar <- LFQServiceData::ionstar$Pep()
#' istar$data <- istar$data %>% dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' filterPep <- LFQService:::filter_proteins_by_peptide_count( istar$data ,  istar$config )
#' tmp <- filter_difference(istar$data, filterPep$data, istar$config)
#' stopifnot(nrow(istar$data )  - nrow(filterPep$data) == nrow(tmp))
#' tmp <- filter_difference(filterPep$data, istar$data , istar$config)
#' stopifnot(nrow(istar$data )  - nrow(filterPep$data) == nrow(tmp))
#'
filter_difference <- function(x, y, config){
  if (nrow(y) > nrow(x)) {
    dplyr::anti_join(y, x, by = config$table$idVars())
  }else{
    dplyr::anti_join(x, y, by = config$table$idVars())
  }
}
