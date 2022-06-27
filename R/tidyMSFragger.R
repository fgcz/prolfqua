#' Methods for reading Fragpipe outputs
#'
#' Convert FragPipe outputs into tidy tables. For more details see functions listed in the see also section.
#'
#' @family FragPipe
#' @name FragPipe
NULL

#'
#' read FragPipe generated MSstats formatted csv files.
#'
#' sanitize entries in the Bioreplicate and Condition columns
#'
#' @family FragPipe
#' @export
#' @param file MSstats formatted file
#' @keywords internal
tidy_FragPipe_MSstats_csv <- function(file){
  inputFile <- readr::read_csv(unz(file, filename = "MSstats.csv"))
  inputFile$BioReplicate <- paste("br", inputFile$BioReplicate, sep = "")
  inputFile$Condition <- make.names(inputFile$Condition)
  inputFile$pep <- 0
  return(inputFile)
}


#' read combined peptide for FragPipe 16 or newer
#' @export
#' @param combprot combinded_peptide.txt data frame read with read.csv.
#' @keywords internal
#' @family FragPipe
tidy_FragPipe_combined_peptides <- function(
  combprot,
  as_list = FALSE
) {
  spcnames = c("spectral.count")
  intnames = c("intensity")
  maxlfqnames = c("maxlfq.intensity")
  protIDcol = "peptide.sequence"

  if (is.character(combprot) && file.exists(combprot)) {
    Cprotein <- as_tibble(read.csv(combprot,
                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE))

  } else if ("tbl_df" %in% class(combprot)) {
    Cprotein <- combprot
  } else {
    stop(class(combprot), " not supported.")
  }
  cnam <- colnames(Cprotein)
  ### start processing
  cnam <- tolower(cnam)
  colnames(Cprotein) <- cnam
  cnam <- cnam[1:which(cnam == "mapped.proteins")]

  message("annotation columns : ", paste(cnam, collapse = "\n"))
  annot <- Cprotein |> dplyr::select(all_of(cnam))
  colnames(Cprotein)

  extractDataLong <- function(Cprotein, what = "total.intensity", butNot = NULL){
    cols <- colnames(Cprotein)
    cols <- setdiff( grep(paste0(what,"$"), cols, value = TRUE) , if (is.null(butNot)) {NULL} else { grep(butNot, cols, value = TRUE) })
    gg <- Cprotein |> dplyr::select( all_of(protIDcol), all_of(cols) )

    gg <- gg |> tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg |> dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", .data$raw.file))
    gg
  }

  res <- vector( mode = "list", length = length(c(intnames, spcnames)))
  names(res)  <- c(intnames, spcnames)

  for (i in seq_along(c(intnames, spcnames))) {
    message("DD: ", c(intnames, spcnames)[i] )
    res[[c(intnames, spcnames)[i]]] <- extractDataLong(Cprotein, what = c(intnames, spcnames)[i], butNot = "maxlfq" )
  }

  if (sum(grepl(".maxlfq.", colnames(Cprotein))) > 0) {
    res_maxlfq <- vector( mode = "list", length(maxlfqnames))
    names(res_maxlfq)  <- maxlfqnames
    for (i in seq_along(maxlfqnames)) {
      message("DD: ", maxlfqnames[i] )
      res_maxlfq[[maxlfqnames[i] ]] <-  extractDataLong(Cprotein, what = maxlfqnames[i], butNot = NULL )
    }
    res <- c(res, res_maxlfq)
  }

  if (as_list) {
    return( res )
  }

  merged <- Reduce( inner_join , res )
  merged <- inner_join( annot, merged )
  return( merged )
}



#' FragPipe read FragPipe combined protein files up to Version 15
#'
#' @export
#' @param combprot path to combined_protein.tsv file
#' @param intnames intensity column prefix
#' @param protIDcol default protein.group
#' @param subgroup default subgroup
#'
#' @keywords internal
#'
#' @family FragPipe
tidy_FragPipe_combined_protein_deprec <- function(
  combprot, intnames = c("total.intensity",
                         "unique.intensity",
                         "razor.intensity",

                         "total.ion.count",
                         "unique.ion.count",
                         "razor.ion.count",

                         "total.spectral.count",
                         "unique.spectral.count",
                         "razor.spectral.count"),
  protIDcol = "protein.group",
  subgroup = "subgroup",
  as_list = FALSE) {
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
  cnam <- cnam[1:which(cnam == "summarized.razor.spectral.count")]
  message("annotation columns : ", paste(cnam, collapse = "\n"))

  annot <- Cprotein |> dplyr::select(all_of(cnam))

  extractDataLong <- function(Cprotein, what = "total.intensity"){
    gg <- Cprotein |> dplyr::select( protIDcol, subgroup, dplyr::ends_with(what))
    gg <- gg |> tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg |> dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", .data$raw.file))
    gg
  }

  res <- vector( mode = "list", length = length(intnames))
  names(res)  <- intnames

  for (i in seq_along(intnames)) {
    res[[intnames[i]]] <- extractDataLong(Cprotein, what = intnames[i] )
  }
  if (as_list) {
    return( res )
  }

  merged <- Reduce(inner_join, res)
  merged <- inner_join(annot, merged)


  return(merged)
}


#' read combined_protein.tsv file for FragPipe Version 16 or newer
#' @export
#' @param combprot path to combined_protein.tsv file
#' @param as_list return as list
#' @return tidy dataframe or list with df (e.g. total.spectral.count or total.intensity etc).
#' @keywords internal
#' @family FragPipe
tidy_FragPipe_combined_protein <- function(
  combprot,
  as_list = FALSE
) {
  spcnames = c("Total Spectral Count",
               "Unique Spectral Count",
               "Razor Spectral Count")
  intnames = c("Total Intensity",
               "Unique Intensity",
               "Razor Intensity")
  maxlfqnames = c("Maxlfq Total Intensity",
                  "Maxlfq Unique Intensity",
                  "Maxlfq Razor Intensity")

  protIDcol = "Protein"

  if (is.character(combprot) && file.exists(combprot)) {
    Cprotein <- as_tibble(read.csv(combprot,
                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))

  } else if ("tbl_df" %in% class(combprot)) {
    Cprotein <- combprot
  } else {
    stop(class(combprot), " not supported.")
  }

  cnam <- gsub("Total Razor ", "Total ",
               gsub("Unique Razor ","Unique ",
                    gsub(" Intensity$"," Razor Intensity",
                         gsub(" Spectral Count$"," Razor Spectral Count",colnames(Cprotein))
                    )
               )
  )

  ### start processing
  cnam <- cnam
  colnames(Cprotein) <- cnam
  cnam <- cnam[1:which(cnam == "Combined Total Spectral Count")]

  message("annotation columns : ", paste(cnam, collapse = "\n"))
  annot <- Cprotein |> dplyr::select(all_of(cnam))
  colnames(Cprotein)

  extractDataLong <- function(Cprotein, what = "Total Intensity", butNot = NULL){
    cols <- colnames(Cprotein)
    cols <- setdiff( grep(paste0(what,"$"), cols, value = TRUE) , if (is.null(butNot)) {NULL} else { grep(butNot, cols, value = TRUE) })
    gg <- Cprotein |> dplyr::select( all_of(protIDcol), all_of(cols) )

    gg <- gg |> tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg |> dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", .data$raw.file))
    gg
  }

  res <- vector( mode = "list", length = length(c(intnames, spcnames)))
  names(res)  <- c(intnames, spcnames)

  for (i in seq_along(c(intnames, spcnames))) {
    message("DD: ", c(intnames, spcnames)[i] )
    res[[c(intnames, spcnames)[i]]] <- extractDataLong(Cprotein, what = c(intnames, spcnames)[i], butNot = "maxlfq" )
  }

  if (sum(grepl(".maxlfq.", colnames(Cprotein))) > 0) {
    res_maxlfq <- vector( mode = "list", length(maxlfqnames))
    names(res_maxlfq)  <- maxlfqnames
    for (i in seq_along(maxlfqnames)) {
      message("DD: ", maxlfqnames[i] )
      res_maxlfq[[maxlfqnames[i] ]] <-  extractDataLong(Cprotein, what = maxlfqnames[i], butNot = NULL )
    }
    res <- c(res, res_maxlfq)
  }

  if (as_list) {
    return( res )
  }

  merged <- Reduce( inner_join , res )
  merged <- inner_join( annot, merged )
  colnames(merged) <- tolower(make.names(colnames(merged)))
  return( merged )
}


