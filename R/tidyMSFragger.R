#' read MSFragger generated MSstats formatted csv files.
#' @export
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

#' read MSFragger combined protein file
#' @export
#' @param
#' @examples
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
tidy_MSFragger_combined_protein <- function(combprot) {
  if (is.character(combprot) && file.exists(combprot)) {
    Cprotein <- as_tibble(read.csv(combprot,
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE))

  } else if("tbl_df" %in% class(combprot)) {
    Cprotein <- combprot
  }

  ### start processing
  colnames(Cprotein) <- tolower(colnames(Cprotein))
  annot <- Cprotein %>% dplyr::select(colnames(Cprotein)[1:14])

  extractDataLong <- function(Cprotein, what = "total.intensity"){
    gg <- Cprotein %>% dplyr::select( protein.group, subgroup, ends_with(what))
    gg <- gg %>% tidyr::pivot_longer(cols = ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg %>% dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", raw.file))
    gg
  }
  intnames <- c("total.intensity",
                "unique.intensity",
                "razor.intensity",
                "total.ion.count",
                "unique.ion.count",
                "razor.ion.count",
                "total.spectral.count",
                "unique.spectral.count",
                "razor.spectral.count")

  res <- vector( mode = "list", length = length(intnames))
  names(res)  <- intnames

  for (i in 1:length(intnames)) {
    res[[intnames[i]]] <- extractDataLong(Cprotein, what = intnames[i] )
  }

  merged <- Reduce(inner_join, res)
  merged <- inner_join(annot, merged)
  return(merged)
}
