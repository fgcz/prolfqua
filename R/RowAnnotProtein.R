
.annotProtein <- function(
    startdata,
    Accession,
    revpattern = "^REV_",
    contpattern = "zz"){

  GRP2 <- list()
  distinctprotid <- startdata |> select(pID = !!sym(Accession)) |> distinct()
  distinctprotid <- distinctprotid |> mutate(
    proteinAnnot = case_when(grepl(revpattern,pID) ~ "REV",
                             grepl(contpattern,pID) ~ "CON",
                             TRUE ~ "FW"))
  GRP2$percentOfContaminants <-  round(mean(distinctprotid$proteinAnnot == "CON") * 100, digits = 2)
  GRP2$percentOfFalsePositives <- round(mean(distinctprotid$proteinAnnot == "REV") * 100, digits = 2)
  GRP2$totalNrOfProteins <- sum( table( distinctprotid$proteinAnnot ) )
  GRP2$NrOfProteinsNoDecoys <- sum(distinctprotid$proteinAnnot == "FW")
  return(list(stats = GRP2, distinctprotid = distinctprotid))
}

# RowAnnotProtein ----
#' Decorates LFQData with a row annotation and some protein specific functions.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#'
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#' pannot <- RowAnnotProtein$new( lfqdata )
#'
#' pannot$annotateREV()
#' pannot$annotateCON()
#'
#' pannot$nr_clean()
#' dd <- pannot$clean()
#' tmp <- lfqdata$get_subset(dd)
#' tmp$complete_cases()
#'
RowAnnotProtein <-
  R6::R6Class("RowAnnotProtein",
              public = list(
                #' @field row_annot data.frame containing further information
                row_annot = NULL,
                #' @field pID column with protein ids
                pID = character(),
                #' @description initialize
                #' @param lfqdata data frame from \code{\link{setup_analysis}}
                #' @param row_annot data frame with row annotation. Must have columns matching \code{config$table$hkeysDepth()}
                initialize = function(lfqdata, row_annot){
                  stopifnot(lfqdata$config$table$hierarchyDepth == 1)
                  self$pID = lfqdata$config$table$hkeysDepth()
                  if (!missing(row_annot)) {
                    row_annot <- dplyr::filter(row_annot, !!sym(self$pID) %in% lfqdata$data[[self$pID]] )
                    stopifnot(self$pID %in% colnames(row_annot))
                    self$row_annot <- row_annot
                  } else {
                    self$row_annot <- distinct(select(lfqdata$data,self$pID))
                  }
                },
                #' @description
                #' annotate rev sequences
                #' @param pattern default "REV_"
                annotateREV = function(pattern = "REV_") {
                  self$row_annot <- self$row_annot |> mutate(
                    REV = case_when(grepl(pattern, !!sym(self$pID), ignore.case = TRUE) ~ TRUE,
                                    TRUE ~ FALSE))

                  return(sum(self$row_annot$REV))
                },
                #' @description
                #' annotate contaminants
                #' @param pattern default "^zz|^CON"
                annotateCON = function(pattern = "^zz|^CON") {
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

