#' Compute correlation matrix
#' @param dataX data.frame with transition intensities per peptide
#' @export
#' @examples
#' data(correlatedPeptideList)
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
#' data(correlatedPeptideList)
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

#' Check if it is a mixed model
#' It checks if the patter (something | something) is present
#' @export
#' @examples
#' model <- "intensity ~ test + (test|test)"
#' is_mixed_model(model)
#' model <- "intensity ~ test"
#' stopifnot(is_mixed_model(model) == FALSE)
is_mixed_model <- function(model){
  res <- grepl("\\(.+\\|.+\\)", model)
  return(res)
}
