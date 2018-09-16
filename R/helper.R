#' id lables for piw
#'
#' @export
#'
getIDLabels <- function(){
  .idlables <- c("Peptide.Sequence", "Protein.Name", "Precursor.Charge", "Product.Charge", "Fragment.Ion", "Isotope.Label")
}
#' piwots light
#' @param data data.frame containing columns : Peptide.Sequence ...
#'@export
piwotPiw <- function(data){
  lightPiw <- reshape2::dcast(data,
                              Peptide.Sequence + Protein.Name + Precursor.Charge +
                                Product.Charge + Fragment.Ion + Isotope.Label
                              ~ Replicate.Name,
                              value.var = "Area")
  idlables <- c("Peptide.Sequence", "Protein.Name", "Precursor.Charge", "Product.Charge", "Fragment.Ion", "Isotope.Label")
  return(lightPiw)
}

#' extracts Intensity columns
#' @param data data.frame returned by piwotPiw
#' @export
getIntensities <- function(data ){
  lightInt <- data[,7:ncol(data)]
  dum <- apply(data[ ,1:6],1, paste,collapse="_")
  rownames(lightInt) <- dum
  return(lightInt)
}


#' Make nice plots of transitions or peptides with correlations
#' @param dataX data.frame
#' @param main some name to plot
#' @param log log transform y axes
#' @param ylab label for y axes
#' @export
#' @examples
#' library(SRMService)
#' library(tidyverse)
#' library(quantable)
#' data(correlatedPeptideList)
#' plotNicely(correlatedPeptideList[[1]])
#'
plotNicely <- function(dataX, main="", log="", ylab="log(intensity)"){
  mat <- matrix(c(1,1,1,1,0,2,2,3), byrow=T, ncol=4)
  graphics::layout(mat, widths=c(2,1,1,1), heights=c(2,1))
  dataXt <- t(dataX)
  graphics::matplot(dataXt,type="l", main=main,lwd=1,lty=1, ylab="log(intensity)",las=2, xaxt = "n", log=log)
  graphics::axis(1, at = 1:nrow(dataXt), labels = rownames(dataXt), cex.axis = 0.7, las=2)
  graphics::legend("bottomleft", legend=rownames(dataX),col=1:5,lty=1 )
  nrow(dataX)
  if(nrow(dataX)>1){
    ordt <- (dataX)[order(apply(dataX,1,mean)),]
    dd <- transitionCorrelationsJack(ordt)
    imageWithLabelsNoLayout(dd,col=getBlueWhiteRed(),zlim=c(-1,1), textB=2)
    imageColorscale(dd,col=getBlueWhiteRed(), zlim=c(-1,1))
    invisible(dd)
  }
}

#' Compute correlation matrix
#' @param dataX data.frame with transition intensities per peptide
#' @export
#' @examples
#' data(correlatedPeptideList)
#' transitionCorrelations(correlatedPeptideList[[1]])
transitionCorrelations <- function(dataX){
  if(nrow(dataX) > 1){
    ordt <- (dataX)[order(apply(dataX,1,mean)),]
    xpep <- t(ordt)
    dd <- stats::cor(t(ordt),use="pairwise.complete.obs", method = "spearman")
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
transitionCorrelationsJack <- function(dataX, distmethod = function(x){cor(x, use="pairwise.complete.obs", method="pearson")}){
  if(nrow(dataX) > 1){
    ordt <- (dataX)[order(apply(dataX,1,mean)),]
    xpep <- t(ordt)
    quantable::jackknifeMatrix(xpep, distmethod)
  }else{
    message("Could not compute correlation, nr rows : " , nrow(dataX) )
  }
}









