% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyMS_R6Model.R
\name{LR_test}
\alias{LR_test}
\title{Likelihood ratio test}
\usage{
LR_test(
  modelProteinF,
  modelName,
  modelProteinF_Int,
  modelName_Int,
  subject_Id = "protein_Id",
  path = NULL
)
}
\arguments{
\item{modelProteinF}{table with models (see build model)}

\item{modelName}{name of model}

\item{modelProteinF_Int}{reduced model}

\item{modelName_Int}{name of reduced model}

\item{subject_Id}{subject id typically Assession or protein_Id}

\item{path}{default NULL, set to a directory if you need to write diagnostic plots.}
}
\description{
Likelihood ratio test
}
\examples{
data_2Factor <- prolfqua::sim_lfq_data_2Factor_config(
 Nprot = 200,
 with_missing = TRUE,
 weight_missing = 2)

pMerged <- LFQData$new(data_2Factor$data, data_2Factor$config)

pMerged$config$table$get_response()
pMerged$factors()

formula_condition_and_Batches <-
  prolfqua::strategy_lm("abundance ~ Treatment + Background")
modCB <- prolfqua::build_model(
  pMerged$data,
  formula_condition_and_Batches,
  subject_Id = pMerged$config$table$hierarchy_keys() )

formula_condition <-
  prolfqua::strategy_lm("abundance ~ Treatment")
modC <- prolfqua::build_model(
  pMerged$data,
  formula_condition,
  subject_Id = pMerged$config$table$hierarchy_keys() )

tmp <- LR_test(modCB$modelDF, "modCB", modC$modelDF, "modB")
hist(tmp$likelihood_ratio_test.pValue)

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
\code{\link{pivot_model_contrasts_2_Wide}()},
\code{\link{plot_lmer_peptide_predictions}()},
\code{\link{sim_build_models_lm}()},
\code{\link{sim_build_models_lmer}()},
\code{\link{sim_make_model_lm}()},
\code{\link{sim_make_model_lmer}()},
\code{\link{strategy_lmer}()},
\code{\link{summary_ROPECA_median_p.scaled}()}
}
\concept{modelling}
