% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6_TransitionCorrelations.R
\name{impute_correlationBased}
\alias{impute_correlationBased}
\title{Imputate missing peptide intensities to maximize peptide correlationion}
\usage{
impute_correlationBased(x, config)
}
\arguments{
\item{x}{data}

\item{config}{configuration}
}
\description{
Assumes the peptide intensities are correlation assumption
}
\examples{



bb <- prolfqua_data('data_ionstar')$normalized()
bb$config <- old2new(bb$config)
config <- bb$config$clone(deep=TRUE)
data <- bb$data
mean(is.na(data$peptide.intensity))
dataI <- impute_correlationBased(data, config)
dim(dataI)
stopifnot(dim(dataI) == c(dim(data)+c(0,1)))
stopifnot(mean(is.na(dataI$srm_ImputedIntensity)) <= mean(is.na(data$peptide.intensity)))

}
\keyword{internal}
