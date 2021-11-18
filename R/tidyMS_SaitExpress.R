
#' Network visualization.
#' look at https://www.jessesadler.com/post/network-analysis-with-r/
#' https://fgcz-intranet.uzh.ch/tiki-index.php?page=WG_APMSnProximityLabeling
#' @examples
#' library(tidyverse)
#'
runSaint <- function(si,
                     filedir = getwd(),
                     spc = TRUE,
                     CLEANUP = TRUE){
  stopifnot(names(si) == c("inter","prey","bait"))
  paths <- character(3)
  for (i in 1:length(si)) {
    filen <- file.path(filedir , paste0(names(si)[i],".txt"))
    paths[i] <- filen
    message(filen)
    readr::write_tsv(si[[i]], file = filen , col_names = FALSE)
  }

  pkg <- find.package("MSProteinInteraction4R")
  if (spc) {
    exeS2 <- "SaintExpress\\bin\\Windows64\\SAINTexpress-spc.exe"
  } else {
    exeS2 <- "SaintExpress\\bin\\Windows64\\SAINTexpress-int.exe"
  }
  exeT <- file.path(pkg, exeS2)

  out <- system2(exeT,
                 args = paths,
                 stdout = TRUE,
                 stderr = TRUE,
                 wait = TRUE,
                 minimized = TRUE)
  Sys.sleep(2) # needed otherwise the list.txt file can't be deleted
  listFile <- file.path(getwd(),"list.txt")
  res <- read.csv(file = listFile, sep = "\t")
  if (CLEANUP) {
    if (!file.remove(listFile)) {
      warning( "can't remove " , listFile)
    }
    file.remove(paths)
  }
  res <- list(listFile = data.frame(listFile = listFile),  list = res, out = data.frame(out = out))
  return(res)
}


#' Convert tidy table with protein quants into saintExpress compatible inputs
#' @export
protein_2localSaint <- function(xx,
                                quantcolumn = "mq.protein.intensity",
                                proteinID = "protein_Id",
                                geneNames  = "protein_Id",
                                proteinLength = "protein.length",
                                IP_name = "raw.file",
                                baitCol = "bait",
                                CorTCol = "CorT"
){
  reqcolumns <- c(quantcolumn,proteinID,geneNames,proteinLength,IP_name,baitCol,CorTCol)
  if ( !all(reqcolumns %in% colnames(xx)) ) {

    stop("columns not found ", paste0(reqcolumns[which(!reqcolumns %in% colnames(xx))]))
  }
  res <- list()
  bait <- xx |> dplyr::select(!!!syms(c(IP_name,baitCol,CorTCol)))
  bait <- dplyr::distinct(bait)
  res$bait <- bait
  prey <- xx |> dplyr::select(!!!syms(c(proteinID,
                                        proteinLength ,geneNames)))
  prey <- distinct(prey)
  res$prey <- prey
  inter <- xx |>
    dplyr::select(!!!syms(c(IP_name,
                            baitCol,
                            proteinID ,
                            quantcolumn))) |>
    filter(!!sym(quantcolumn) > 0)
  res$inter <- inter
  res <- res[c("inter","prey","bait")]
  return(res)
}




#' add protein lengths from fasta file to data frame (id_col - protein id column.)
#' @export
addProteinLengths <- function(intdata, fasta_file , id_col = "protein_Id" ){
  fasta <- prozor::readPeptideFasta(file = fasta_file)
  plengths <- data.frame(id = names(fasta) , protein.length = sapply(fasta, stringr::str_length))
  byx <- "id"
  names(byx) <- id_col
  intdata <- dplyr::left_join(intdata, plengths , by = byx)
  intdata$protein.length[is.na(intdata$protein.length)] <- as.integer(mean(intdata$protein.length, na.rm = TRUE))
  return(intdata)
}
