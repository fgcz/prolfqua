# ProteinAnnotation ----
#' Decorates LFQData with a row annotation and some protein specific functions.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <-prolfqua::sim_lfq_data_peptide_config()
#' lfqdata <- LFQData$new(istar$data, istar$config)
#' pannot <- ProteinAnnotation$new( lfqdata )
#' pannot$annotate_decoys()
#' pannot$annotate_contaminants()
#' dd <- pannot$clean()
#' tmp <- lfqdata$get_subset(dd)
#' pannot$row_annot
#'
ProteinAnnotation <-
  R6::R6Class("ProteinAnnotation",
              public = list(
                #' @field row_annot data.frame containing further information
                row_annot = NULL,
                #' @field pID column with protein ids
                pID = character(),
                #' @field description name of column containing descriptions
                description = "description",
                #' @field ids vector with columns containing addition IDs
                ids = character(),
                #' @field nr_children name of columns with the number of peptides
                nr_children = character(),
                #' @description initialize
                #' @param lfqdata data frame from \code{\link{setup_analysis}}
                #' @param row_annot data frame with row annotation. Must have columns matching \code{config$table$hierarchy_keys_depth()}
                #' @param description name of column with description
                #' @param ids names of columns with additional ID's
                #' @param nr_peptides additional peptides
                initialize = function(lfqdata,
                                      row_annot,
                                      description = NULL,
                                      ids = NULL,
                                      nr_children = "nr_peptides"){
                  self$pID = lfqdata$config$table$hierarchy_keys_depth()[1]
                  self$nr_children = nr_children
                  if (!missing(row_annot)) {
                    stopifnot(self$pID %in% colnames(row_annot))
                    row_annot <- dplyr::filter(row_annot, !!sym(self$pID) %in% lfqdata$data[[self$pID]] )
                    self$row_annot <- row_annot
                  } else {
                    self$row_annot <- distinct(select(lfqdata$data, self$pID))
                  }
                  if (!self$nr_children %in% colnames(row_annot) ) {
                    self$row_annot <- inner_join(
                      self$row_annot,
                      nr_obs_experiment(lfqdata$data, lfqdata$config, name_nr_child = self$nr_children),
                      by = self$pID)
                  }
                },
                #' @description
                #' annotate rev sequences
                #' @param pattern default "REV_"
                annotate_decoys = function(pattern = "REV_") {
                  self$row_annot <- self$row_annot |> mutate(
                    REV = case_when(grepl(pattern, !!sym(self$pID), ignore.case = TRUE) ~ TRUE,
                                    TRUE ~ FALSE))

                  return(sum(self$row_annot$REV))
                },
                #' @description
                #' annotate contaminants
                #' @param pattern default "^zz|^CON"
                annotate_contaminants = function(pattern = "^zz|^CON") {
                  self$row_annot <- self$row_annot |> mutate(
                    CON = case_when(grepl(pattern, !!sym(self$pID), ignore.case = TRUE) ~ TRUE,
                                    TRUE ~ FALSE))
                  return(sum(self$row_annot$CON))
                },
                #' @description get number of neither contaminants nor decoys
                #' @param contaminants remove contaminants
                #' @param decoys remove decoys
                #' return number of cleans
                nr_clean = function(contaminants = TRUE, decoys = TRUE){

                  if (decoys && !("REV" %in% colnames(self$row_annot)) ) { stop("annotate REV") }
                  if (contaminants & !("CON" %in% colnames(self$row_annot)) ) { stop("annotate CON") }

                  res <- if (decoys && contaminants) {
                    sum(!self$row_annot$REV & !self$row_annot$CON)
                  } else if (contaminants) {
                    sum(!self$row_annot$CON)
                  } else if (decoys) {
                    sum(!self$row_annot$REV)
                  } else {
                    nrow(self$row_annot)
                  }
                  return(res)
                },
                #' @description remove REV and CON sequences
                #' @param contaminants remove contaminants
                #' @param decoys remove decoys
                #'
                clean = function(contaminants = TRUE, decoys = TRUE){
                  if (contaminants && !("REV" %in% colnames(self$row_annot)) ) { stop("annotate REV") }
                  if (decoys && !("CON" %in% colnames(self$row_annot)) ) { stop("annotate CON") }
                  res <- if (decoys && contaminants) {
                    filter(self$row_annot , !self$row_annot$REV & !self$row_annot$CON )
                  } else if (contaminants) {
                    filter(self$row_annot , !self$row_annot$CON)
                  } else if (decoys) {
                    filter(self$row_annot , !self$row_annot$REV )
                  } else {
                    self$row_annot
                  }
                  return(res)
                }

              )
  )

