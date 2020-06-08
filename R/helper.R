#' Compute correlation matrix
#' @param dataX data.frame with transition intensities per peptide
#' @export
#' @examples
#' correlatedPeptideList <- LFQServiceData::correlatedPeptideList
#' transitionCorrelations(correlatedPeptideList[[1]])
transitionCorrelations <- function(dataX) {
  if (nrow(dataX) > 1) {
    ordt <- (dataX)[order(apply(dataX, 1, mean)), ]
    dd <- stats::cor(t(ordt), use = "pairwise.complete.obs", method = "spearman")
    return(dd)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
  }
}

#' Compute correlation matrix with jack
#' @param dataX data.frame with transition intensities per peptide
#' @export
#' @importFrom stats cor
#' @examples
#' correlatedPeptideList <- LFQServiceData::correlatedPeptideList
#' transitionCorrelationsJack(correlatedPeptideList[[1]])
transitionCorrelationsJack <- function(dataX,
                                       distmethod =
                                         function(x){cor(x, use = "pairwise.complete.obs", method = "pearson")}) {
  if (nrow(dataX) > 1) {
    ordt <- (dataX)[order(apply(dataX, 1, mean)), ]
    xpep <- t(ordt)
    LFQService::jackknifeMatrix(xpep, distmethod)
  }else{
    message("Could not compute correlation, nr rows : ", nrow(dataX))
  }
}

