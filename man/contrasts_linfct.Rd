% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6_Modelling.R
\name{contrasts_linfct}
\alias{contrasts_linfct}
\title{compute group averages}
\usage{
contrasts_linfct(
  models,
  linfct,
  subject_Id = "protein_Id",
  contrastfun = prolfqua::my_contest
)
}
\description{
used in p2621, p2109
}
\examples{
modelSummary_A <- sim_build_models_lm()
m <- get_complete_model_fit(modelSummary_A$modelDF)

factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
factor_levelContrasts <- contrasts_linfct( m,
        factor_contrasts,
        subject_Id = "protein_Id",
        contrastfun = prolfqua::my_contrast_V2)


}
\keyword{internal}
