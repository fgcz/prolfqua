#' read MSFragger generated MSstats formatted csv files.
#' @export
#' @rdname MSFragger
#' @param file MSstats formatted file
#' @family MSFragger
#' @keywords internal
#'
tidy_MSFragger_MSstats_csv <- function(file){
  inputFile <- readr::read_csv(unz(file, filename = "MSstats.csv"))
  inputFile$BioReplicate <- paste("br", inputFile$BioReplicate, sep = "")
  inputFile$Condition <- make.names(inputFile$Condition)
  inputFile$pep <- 0
  return(inputFile)
}




#' read MSFragger combined protein file up to Version 15
#' @export
#' @rdname MSFragger
#' @param combprot path to combined_protein.tsv file
#' @param intnames intensity column prefix
#' @family MSFragger
#' @keywords internal
#' @examples
#'
#' if(FALSE){
#'   unzip(inputMQfile, list = TRUE)$Name
#'   protein <- as_tibble(read.csv(unz(inputMQfile,"IonstarWithMSFragger/combined_protein.tsv"),
#'                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE))
#'   tidy_MSFragger_combined_protein(protein)
#' }
#'
tidy_MSFragger_combined_protein <- function(combprot, intnames = c("total.intensity",
                                                                   "unique.intensity",
                                                                   "razor.intensity",

                                                                   "total.ion.count",
                                                                   "unique.ion.count",
                                                                   "razor.ion.count",

                                                                   "total.spectral.count",
                                                                   "unique.spectral.count",
                                                                   "razor.spectral.count"),
                                            protIDcol = "protein.group", subgroup = "subgroup") {
  if (is.character(combprot) && file.exists(combprot)) {
    Cprotein <- as_tibble(read.csv(combprot,
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE))

  } else if ("tbl_df" %in% class(combprot)) {
    Cprotein <- combprot
  } else {
    stop(class(combprot), " not supported.")
  }

  ### start processing
  colnames(Cprotein) <- tolower(colnames(Cprotein))
  cnam <- tolower(colnames(Cprotein))
  cnam <- cnam[1:which(cnam == "combined.total.spectral.count")]
  message("annotation columns : ", paste(cnam, collapse = "\n"))

  annot <- Cprotein %>% dplyr::select(all_of(cnam))

  extractDataLong <- function(Cprotein, what = "total.intensity"){
    gg <- Cprotein %>% dplyr::select( protIDcol, subgroup, dplyr::ends_with(what))
    gg <- gg %>% tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg %>% dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", .data$raw.file))
    gg
  }

  res <- vector( mode = "list", length = length(intnames))
  names(res)  <- intnames

  for (i in seq_along(intnames)) {
    res[[intnames[i]]] <- extractDataLong(Cprotein, what = intnames[i] )
  }
  return(res)
  merged <- Reduce(inner_join, res)
  merged <- inner_join(annot, merged)
  return(merged)
}


#' read MSFragger combined protein file
#' @export
#' @rdname MSFragger
#' @param combprot path to combined_protein.tsv file
#' @param as_list return as list
#' @return tidy dataframe or list with df (e.g. total.spectral.count or total.intensity etc).
#' @family MSFragger
#' @keywords internal
#' @examples
#'
#' if(FALSE){
#'   unzip(inputMQfile, list = TRUE)$Name
#'   protein <- as_tibble(read.csv(unz(inputMQfile,"IonstarWithMSFragger/combined_protein.tsv"),
#'                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE))
#'   tidy_MSFragger_combined_protein(protein)
#' }
#'
tidy_MSFragger_combined_protein_V16 <- function(
  combprot,

  as_list = FALSE
) {
  spcnames = c("total.spectral.count",
               "unique.spectral.count",
               "razor.spectral.count")
  intnames = c("total.intensity",
               "unique.intensity",
               "razor.intensity")
  maxlfqnames = c("maxlfq.total.intensity",
                  "maxlfq.unique.intensity",
                  "maxlfq.razor.intensity")
  protIDcol = "protein"

  if (is.character(combprot) && file.exists(combprot)) {
    Cprotein <- as_tibble(read.csv(combprot,
                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE))

  } else if ("tbl_df" %in% class(combprot)) {
    Cprotein <- combprot
  } else {
    stop(class(combprot), " not supported.")
  }

  Cprotein <- protein
  cnam <- gsub("Total.Razor.", "Total.",
               gsub("Unique.Razor.","Unique.",
                    gsub("\\.Intensity$",".Razor.Intensity",
                         gsub("\\.Spectral.Count$",".Razor.Spectral.Count",colnames(Cprotein))
                    )
               )
  )

  ### start processing
  cnam <- tolower(cnam)
  colnames(Cprotein) <- cnam
  cnam <- cnam[1:which(cnam == "combined.total.spectral.count")]

  message("annotation columns : ", paste(cnam, collapse = "\n"))
  annot <- Cprotein %>% dplyr::select(all_of(cnam))
  colnames(Cprotein)

  extractDataLong <- function(Cprotein, what = "total.intensity", butNot = NULL){
    cols <- colnames(Cprotein)
    cols <- setdiff( grep(paste0(what,"$"), cols, value = TRUE) , if (is.null(butNot)) {NULL} else { grep(butNot, cols, value = TRUE) })
    gg <- Cprotein %>% dplyr::select( all_of(protIDcol), all_of(cols) )

    gg <- gg %>% tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg %>% dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", .data$raw.file))
    gg
  }

  res <- vector( mode = "list", length = length(c(intnames, spcnames, maxlfqnames)))
  names(res)  <- c(intnames, spcnames, maxlfqnames)

  for (i in seq_along(c(intnames, spcnames))) {
    message("DD: ", c(intnames, spcnames)[i] )
    res[[c(intnames, spcnames)[i]]] <- extractDataLong(Cprotein, what = c(intnames, spcnames)[i], butNot = "maxlfq" )
  }

  for (i in seq_along(maxlfqnames)) {
    message("DD: ", maxlfqnames[i] )
    res[[maxlfqnames[i] ]] <-  extractDataLong(Cprotein, what = maxlfqnames[i], butNot = NULL )
  }

  if (as_list) {
    return(res)
  }
  merged <- Reduce(inner_join, res)
  merged <- inner_join(annot, merged)
  return(merged)
}
