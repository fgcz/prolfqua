% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6_TransitionCorrelations.R
\name{normalize_log2_robscale}
\alias{normalize_log2_robscale}
\title{normalize data by log2 and robust scaling}
\usage{
normalize_log2_robscale(pdata, config)
}
\arguments{
\item{pdata}{data.frame}

\item{config}{AnalysisConfiguration}
}
\value{
list with data.frame (data) and updated config (config)
}
\description{
normalize data by log2 and robust scaling
}
\examples{

istar <- prolfqua_data('data_ionstar')$filtered()
istar$config <- old2new(istar$config)
istar_data <- istar$data |> dplyr::filter(protein_Id \%in\% sample(protein_Id, 100))
xx <- normalize_log2_robscale(istar_data, istar$config)
names(xx)
xx$config$table$workIntensity

}
\seealso{
Other preprocessing: 
\code{\link{INTERNAL_FUNCTIONS_BY_FAMILY}},
\code{\link{apply_to_response_matrix}()},
\code{\link{filter_proteins_by_peptide_count}()},
\code{\link{get_robscales}()},
\code{\link{robust_scale}()},
\code{\link{scale_with_subset}()},
\code{\link{scale_with_subset_by_factors}()}
}
\concept{preprocessing}
\keyword{internal}
