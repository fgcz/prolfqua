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
      dplyr::mutate(!!paste0("factor_",factor) := purrr::map(!!sym(paste0("lmer_",modelName )), ~.contrast_tukey_multcomp(.,factor=factor)))
  }
  dd <- modelProteinF %>% dplyr::select(subject_Id, starts_with("factor_"))
  contrasts <- dd %>% tidyr::gather("factor", "contrasts", - c(!!!(syms( subject_Id)))) %>%
    tidyr::unnest() %>% arrange(!!!syms(subject_Id)) %>%
    dplyr::select(-rhs)

  return(contrasts)
}


#' compute all contrasts from non interaction model automatically.
#'
#' used p2109
#' @export
#' @import tidyverse
#' @import magrittr
#' @examples
#' library(tidyverse)
#' D <- LFQService::resultsV12954
#' formula_randomPeptide <- make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)")
#' modelName <- "f_Condition_r_peptide"
#' pepConfig <- D$config_pepIntensityNormalized
#' modellingResult <- workflow_model_analyse( D$pepIntensityNormalized,
#' formula_randomPeptide,
#' subject_Id = pepConfig$table$hkeysLevel(),
#' modelName)
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
  contrasts <- dplyr::inner_join(modelProteinF, contrasts )

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


.setLargeQValuesToNA <- function(data,
                                 QValueColumn,
                                 intensityOld,
                                 thresholdQValue = 0.05,
                                 intensityNew = "IntensitiesWithNA"){

  thresholdF <- function(x,y, threshold = 0.05){ ifelse(x < threshold, y, NA)}
  data <- data %>%
    dplyr::mutate(!!intensityNew := thresholdF(!!!syms(c(QValueColumn ,intensityOld )), threshold = thresholdQValue))
  return(data)
}


#' Perform anova analysis
#' @export
#' @importFrom glue glue
#' @importFrom purrr map
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' library(glue)
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' data <- transform_work_intensity(data, config, log2)
#' x1 <- deprecated_compute_anova_lm(data, config, hierarchy_level= 2, factor_level=1)
#' x2 <- deprecated_compute_anova_lm(data, config, hierarchy_level = 1, factor_level=2)
#' x3 <- deprecated_compute_anova_lm(data, config, hierarchy_level= 2, factor_level=2)
#' (head(x1))
deprecated_compute_anova_lm <- function(data, config, .formula=NULL, hierarchy_level=1, factor_level=1){
  aovmodelfit <- function(x, formula){
    tryCatch(anova(lm(formula , data=x)), error = function(e) return(NULL))
  }

  if(is.null(.formula)){
    formulastr <- paste(config$table$getWorkIntensity(), " ~ ", paste(c(config$table$factorKeys()[1:factor_level],
                                                                        config$table$hierarchyKeys(TRUE)[1],
                                                                        config$table$fileName), collapse=" + "))
    print(formulastr)
    formula <- as.formula(formulastr)
  }
  message("formula :" , deparse(formula))
  groupVars <-config$table$hierarchyKeys()[1:hierarchy_level]

  pepRes <- data %>% dplyr::group_by(!!!syms(groupVars)) %>% tidyr::nest()
  pepRes1 <- pepRes %>% dplyr::mutate(anova = map( data, aovmodelfit, formula))
  head(pepRes1)
  pepRes2 <- pepRes1 %>% dplyr::mutate(broomres = map(anova, broom::tidy))
  pepRes3 <- pepRes2 %>%
    dplyr::select(!!!syms(c(groupVars, "broomres"))) %>%
    tidyr::unnest() %>%
    dplyr::filter(term != "Residuals")

  pVals <- pepRes3 %>% dplyr::select(!!!syms(c(groupVars,"term","p.value")))

  pVals <- pVals %>% dplyr::mutate(term = glue("{term}.p.value"))
  pVals <- pVals %>% tidyr::spread("term", "p.value")
  df <- pepRes3 %>% dplyr::select(!!!syms(c(groupVars,"term","df")))
  df <- df %>% dplyr::mutate(term = glue("{term}.df"))
  df <- df %>% tidyr::spread(term, df)
  statistic <- pepRes3 %>% dplyr::select(!!!syms(c(groupVars,"term","statistic")))
  statistic <- statistic %>% dplyr::mutate(term = glue("{term}.statistic"))
  statistic <- statistic %>% tidyr::spread(term, statistic)
  res <- dplyr::inner_join(dplyr::inner_join(pVals, df, by=groupVars), statistic, by=groupVars)
  return(res)
}


