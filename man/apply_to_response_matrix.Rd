% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6_TransitionCorrelations.R
\name{apply_to_response_matrix}
\alias{apply_to_response_matrix}
\title{Apply function requiring a matrix to tidy table}
\usage{
apply_to_response_matrix(data, config, .func, .funcname = NULL)
}
\arguments{
\item{data}{data.frame}

\item{config}{AnalysisConfiguration}

\item{.func}{function}

\item{.funcname}{name of function (used for creating new column)}
}
\description{
Apply function requiring a matrix to tidy table
}
\examples{


bb <- sim_lfq_data_peptide_config(Nprot = 100)
data <- bb$data
conf <- bb$config
res <- apply_to_response_matrix(data, conf, .func = base::scale)

stopifnot("abundance_base..scale" \%in\% colnames(res))
stopifnot("abundance_base..scale" == conf$table$get_response())
conf <- bb$config$clone(deep=TRUE)
conf$table$workIntensity <- "abundance"
res <- apply_to_response_matrix(data, conf$clone(deep=TRUE), .func = robust_scale)

# Normalize data using the vsn method from bioconductor

if( require("vsn")){
 res <- apply_to_response_matrix(data, conf$clone(deep=TRUE), .func = vsn::justvsn)
}

}
\seealso{
Other preprocessing: 
\code{\link{INTERNAL_FUNCTIONS_BY_FAMILY}},
\code{\link{filter_proteins_by_peptide_count}()},
\code{\link{get_robscales}()},
\code{\link{normalize_log2_robscale}()},
\code{\link{robust_scale}()},
\code{\link{scale_with_subset}()},
\code{\link{scale_with_subset_by_factors}()}
}
\concept{preprocessing}
\keyword{internal}
