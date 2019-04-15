# get all comparisons
.contrast_tukey_multcomp <- function(model, factor){
  if(class(model) == "lm") # fixes issue of mutlcomp not working on factors of class character
  {
    model$model <- as.data.frame(unclass(model$model))
  }

  mcpDef <- multcomp::mcp(Dummy="Tukey")
  names(mcpDef) <- factor
  glt <- multcomp::glht(model, mcpDef)
  sglt <- summary(glt)
  sglt <- broom::tidy(sglt)
  ciglt <- broom::tidy(confint(glt)) %>% dplyr::select(-estimate)
  xx <- dplyr::inner_join(
    sglt,
    ciglt, by=c("lhs","rhs")
  )
  return(xx)
}


# compute contrasts for all factors.
.compute_contrasts_no_interaction <- function( modelProteinF, pepConfig, modelName){
  factors <- pepConfig$table$factorKeys()[1:pepConfig$table$factorLevel]
  subject_Id <- pepConfig$table$hkeysLevel()

  for(factor in factors){
    print(factor)
    modelProteinF <- modelProteinF %>%
      mutate(!!paste0("factor_",factor) := purrr::map(!!sym(paste0("lmer_",modelName )), ~.contrast_tukey_multcomp(.,factor=factor)))
  }
  dd <- modelProteinF %>% dplyr::select(subject_Id, starts_with("factor_"))
  contrasts <- dd %>% gather("factor", "contrasts", - c(!!!(syms( subject_Id)))) %>% unnest() %>% arrange(!!!syms(subject_Id)) %>% dplyr::select(-rhs)

  return(contrasts)
}


#' compute all contrasts from non interaction model automatically.
#'
#' used p2109
#' @export
#' @examples
#' D <- LFQService::resultsV12954
#' formula_randomPeptide <- make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)")
#' modelName <- "f_Condition_r_peptide"
#' pepConfig <- D$config_pepIntensityNormalized
#' #modellingResult <- workflow_interaction_modelling(D$pepIntensityNormalized, D$config_pepIntensityNormalized,formula_randomPeptide, modelName)
#' modellingResult <- workflow_model_analyse( D$pepIntensityNormalized,
#' D$config_pepIntensityNormalized,
#' formula_randomPeptide,
#' modelName,
#' isSingular = lme4::isSingular)
#' results <- deprecated_model_contrasts_no_interaction(modellingResult$modelProteinF,
#' modelName,
#' D$config_pepIntensityNormalized)
#'
#' if( FALSE ){
#'  results$fig$VolcanoPlot
#'  results$fig$histogram_coeff_p.values
#' }
deprecated_model_contrasts_no_interaction <- function(modelProteinF,
                                                    modelName,
                                                    pepConfig
)
{
  warning("Softdeprecated - use explicit linear functions and workflow_model_contrasts function")
  # TODO make chack that model has no contrasts!
  result <- list()

  contrasts <- .compute_contrasts_no_interaction(modelProteinF, pepConfig, modelName )
  modelProteinF<- modelProteinF %>% dplyr::select_at(c( pepConfig$table$hkeysLevel(), "isSingular"))
  contrasts <- inner_join(modelProteinF, contrasts )

  #result$contrasts <- contrasts
  result$fig <- list()
  result$fig$histogram_coeff_p.values_name <- paste0("Contrasts_Auto_histogram_",modelName,".pdf")
  result$fig$histogram_coeff_p.values <- ggplot(data=contrasts,
                                                aes(x = p.value, group=lhs)) + geom_histogram(bins = 20) + facet_wrap(~lhs)


  result$fig$VolcanoPlot_name <- paste0("Contrasts_Auto_Volcano_",modelName,".pdf")
  result$fig$VolcanoPlot <- quantable::multigroupVolcano(contrasts,
                                                         effect = "estimate",
                                                         type = "p.value",
                                                         condition = "lhs",
                                                         label = pepConfig$table$hkeysLevel(),
                                                         xintercept = c(-1, 1),colour = "isSingular")
  result$contrasts <- contrasts
  return(result)
}

