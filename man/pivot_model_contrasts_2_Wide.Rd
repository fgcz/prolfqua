% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6_Modelling.R
\name{pivot_model_contrasts_2_Wide}
\alias{pivot_model_contrasts_2_Wide}
\title{pivot model contrasts matrix to wide format produced by `contrasts_linfct` and ...}
\usage{
pivot_model_contrasts_2_Wide(
  modelWithInteractionsContrasts,
  subject_Id = "protein_Id",
  columns = c("estimate", "p.value", "p.value.adjusted"),
  contrast = "lhs"
)
}
\description{
pivot model contrasts matrix to wide format produced by `contrasts_linfct` and ...
}
\examples{

# this function is used by the contrast classes to implement the to wide method

}
\seealso{
Other modelling: 
\code{\link{Contrasts}},
\code{\link{ContrastsMissing}},
\code{\link{ContrastsModerated}},
\code{\link{ContrastsPlotter}},
\code{\link{ContrastsProDA}},
\code{\link{ContrastsROPECA}},
\code{\link{ContrastsTable}},
\code{\link{INTERNAL_FUNCTIONS_BY_FAMILY}},
\code{\link{LR_test}()},
\code{\link{Model}},
\code{\link{build_model}()},
\code{\link{contrasts_fisher_exact}()},
\code{\link{get_anova_df}()},
\code{\link{get_complete_model_fit}()},
\code{\link{get_p_values_pbeta}()},
\code{\link{isSingular_lm}()},
\code{\link{linfct_all_possible_contrasts}()},
\code{\link{linfct_factors_contrasts}()},
\code{\link{linfct_from_model}()},
\code{\link{linfct_matrix_contrasts}()},
\code{\link{merge_contrasts_results}()},
\code{\link{model_analyse}()},
\code{\link{model_summary}()},
\code{\link{moderated_p_limma}()},
\code{\link{moderated_p_limma_long}()},
\code{\link{my_contest}()},
\code{\link{my_contrast}()},
\code{\link{my_contrast_V1}()},
\code{\link{my_contrast_V2}()},
\code{\link{my_glht}()},
\code{\link{plot_lmer_peptide_predictions}()},
\code{\link{sim_build_models_lm}()},
\code{\link{sim_build_models_lmer}()},
\code{\link{sim_make_model_lm}()},
\code{\link{sim_make_model_lmer}()},
\code{\link{strategy_lmer}()},
\code{\link{summary_ROPECA_median_p.scaled}()}
}
\concept{modelling}
\keyword{internal}
