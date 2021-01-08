#' read MSFragger generated MSstats formatted csv files.
#' @export
#' @param file MSstats formatted file
#'
read_MSFragger_MSstats_csv <- function(file){
  inputFile <- readr::read_csv(unz(file, filename = "MSstats.csv"))
  inputFile$BioReplicate <- paste("br", inputFile$BioReplicate, sep = "")
  inputFile$Condition <- make.names(inputFile$Condition)
  inputFile$pep <- 0
  return(inputFile)
}
