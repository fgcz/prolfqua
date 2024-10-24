% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_LFQ_data.R
\name{sim_lfq_data_protein_config}
\alias{sim_lfq_data_protein_config}
\title{Simulate data, protein, with config}
\usage{
sim_lfq_data_protein_config(
  Nprot = 10,
  with_missing = TRUE,
  weight_missing = 0.2,
  seed = 1234,
  paired = FALSE
)
}
\arguments{
\item{with_missing}{add missing values, default TRUE}

\item{seed}{seed for reproducibility, if NULL no seed is set.}

\item{description}{Nprot number of proteins}
}
\description{
Simulate data, protein, with config
}
\examples{

x <- sim_lfq_data_protein_config()
stopifnot("data.frame" \%in\% class(x$data))
stopifnot("AnalysisConfiguration" \%in\% class(x$config))
x <- sim_lfq_data_protein_config(with_missing = FALSE)

stopifnot(sum(is.na(x$data$abundance)) == 0)
# debug(sim_lfq_data_protein_config)
xp <- sim_lfq_data_protein_config(with_missing = FALSE, paired = TRUE)
stopifnot(length(xp$config$table$factors) == 2)
stopifnot(nrow(xp$data) == nrow(x$data))
}
