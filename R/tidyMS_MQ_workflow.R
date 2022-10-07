# Helper functions -----

#' Keep only those proteins with 2 IDENTIFIED peptides
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @return list with data.frame (data) and name of new column (name)
#' @export
#' @keywords internal
#' @family preprocessing
#' @examples
#'
#'
#'
#' istar <- prolfqua_data('data_ionstar')$Pep()
#' istar$config <- old2new(istar$config)
#' istar_data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' filterPep <- prolfqua::filter_proteins_by_peptide_count( istar_data ,  istar$config )
#'  x <- prolfqua::summarize_hierarchy(filterPep$data , istar$config)
#' stopifnot(x$peptide_Id_n >= istar$config$parameter$min_peptides_protein)
#'
#'
filter_proteins_by_peptide_count <-
  function(pdata,
           config){

    # remove single hit wonders
    tmp <- prolfqua::nr_B_in_A(pdata,config)
    if (!is.null(tmp)) {
      res <- dplyr::filter(tmp$data, !!sym(tmp$name) >= config$parameter$min_peptides_protein )
      name <- tmp$name
    }else{
      res <- pdata
      name <- NULL
    }
    return(list(data = res, name = name))
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
#'
#'
#' istar <- prolfqua_data('data_ionstar')$Pep()
#' istar$config <- old2new(istar$config)
#' istar_data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' filterPep <- prolfqua:::filter_proteins_by_peptide_count( istar_data ,  istar$config )
#' tmp <- filter_difference(istar_data, filterPep$data, istar$config)
#' stopifnot(nrow(istar_data )  - nrow(filterPep$data) == nrow(tmp))
#' tmp <- filter_difference(filterPep$data, istar_data , istar$config)
#' stopifnot(nrow(istar_data )  - nrow(filterPep$data) == nrow(tmp))
#'
filter_difference <- function(x, y, config){
  if (nrow(y) > nrow(x)) {
    dplyr::anti_join(y, x, by = config$table$id_vars())
  }else{
    dplyr::anti_join(x, y, by = config$table$id_vars())
  }
}


