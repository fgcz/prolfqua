# ProteinAnnotation ----
#' Decorates LFQData with a row annotation and some protein specific functions.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' istar$config <- old2new(istar$config)
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' pannot <- ProteinAnnotation$new( lfqdata )
#'
#' pannot$annotate_decoys()
#' pannot$annotate_contaminants()
#'
#' pannot$nr_clean()
#' dd <- pannot$clean()
#' tmp <- lfqdata$get_subset(dd)
#' tmp$complete_cases()
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
                #' @field nr_peptides name of columns with the number of peptides
                nr_peptides = character(),
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
                                      nr_peptides = "nr_peptides"){
                  self$pID = lfqdata$config$table$hierarchy_keys_depth()[1]
                  if (!missing(row_annot)) {
                    stopifnot(self$pID %in% colnames(row_annot))
                    row_annot <- dplyr::filter(row_annot, !!sym(self$pID) %in% lfqdata$data[[self$pID]] )
                    self$row_annot <- row_annot
                  } else {
                    self$row_annot <- distinct(select(lfqdata$data, self$pID))
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
                #' @description
                #' return number of cleans
                nr_clean = function(){
                  if (!("REV" %in% colnames(self$row_annot)) ) { stop("annotate REV") }
                  if (!("CON" %in% colnames(self$row_annot)) ) { stop("annotate CON") }
                  return(sum(!self$row_annot$REV & !self$row_annot$CON))
                },
                #' @description
                #' remove REV and CON sequences
                clean = function(){
                  if (!("REV" %in% colnames(self$row_annot)) ) { stop("annotate REV") }
                  if (!("CON" %in% colnames(self$row_annot)) ) { stop("annotate CON") }
                  return(filter(self$row_annot , !self$row_annot$REV & !self$row_annot$CON) )
                }

              )
  )

