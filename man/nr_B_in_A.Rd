% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6_TransitionCorrelations.R
\name{nr_B_in_A}
\alias{nr_B_in_A}
\title{Compute nr of B per A}
\usage{
nr_B_in_A(pdata, config, merge = TRUE)
}
\arguments{
\item{pdata}{data.frame}

\item{config}{AnalysisConfiguration}
}
\description{
Compute nr of B per A
}
\examples{

bb <- prolfqua_data('data_ionstar')$filtered()
bb$config <- old2new(bb$config)
stopifnot(nrow(bb$data) == 25780)
config <- bb$config$clone(deep=TRUE)
data <- bb$data |> dplyr::select(-all_of("nr_peptide_Id_IN_protein_Id"))
hierarchy <- config$table$hierarchy_keys()
res <- nr_B_in_A(data, config)

res$data |>
  dplyr::select_at(c(config$table$hierarchy_keys_depth(),  res$name)) |>
  dplyr::distinct() |>
  dplyr::pull() |> table()


bb <- prolfqua_data('data_skylineSRM_HL_A')
config <- old2new(bb$config_f())
data <- bb$data
data$Area[data$Area == 0] <- NA
analysis <- setup_analysis(data, config)

resDataStart <- bb$analysis(bb$data, config)


nr_B_in_A(resDataStart, config)
nr_B_in_A(resDataStart, config, merge = FALSE)
config$table$hierarchyDepth <- 2
nr_B_in_A(resDataStart, config, merge = FALSE)

bb <- prolfqua_data('data_IonstarProtein_subsetNorm')
bb$config <- old2new(bb$config$clone(deep=TRUE))
nr_B_in_A(bb$data, bb$config)
#undebug(nr_B_in_A)
}
\keyword{internal}
