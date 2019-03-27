# function for modelling go here.
#' rocs helper function
rocs <- function(data ,response, predictor){
  responseX <- data %>% pull(!!sym(response))
  predictorX <- data %>% pull(!!sym(predictor))
  levels = levels(as.factor(responseX))
  if(length(levels) < 2){
    return(NULL)
  }
  res <- list()
  comparisons <- combn(levels, 2)
  for(i in 1:ncol(comparisons)){
    comp <- comparisons[,i]
    res[[i]] <-  tryCatch(pROC::roc(response = responseX,
                                    predictor = predictorX , levels=comp), error = function(x) NULL)
  }
  return(res)
}

#' Apply roc analysis on main factor on lowest level
#' @export
#' @importFrom purrr map
#' @examples
#' library(tidyverse)
#' library(LFQService)
#' config <- spectronautDIAData250_config$clone(deep=TRUE)
#' config$parameter$min_nr_of_notNA  <- 20
#' data <- spectronautDIAData250_analysis
#' x <- sample(data$protein_Id,2)
#' data <- data %>% dplyr::filter(protein_Id %in% x)
#' res <- compute_roc(data, config)
#' head(res)
#' i <- 2
#'
#' pROC::plot.roc(res$rocs[[i]], print.auc = TRUE, main = paste(res$protein_Id[[i]], "\n",paste(res$rocs[[i]]$levels, collapse = " vs ")))
#' unique(res$protein_Id)
compute_roc <- function(data, config){
  nested <- data %>% group_by(!!sym(config$table$hierarchyKeys()[1]) ,
                              !!sym(config$table$hierarchyKeys(TRUE)[1])) %>% nest()
  nested <- nested %>% mutate(rocs = map(data ,
                                         rocs, response = config$table$factorKeys()[1],
                                         predictor= config$table$getWorkIntensity() ))

  nested <- nested %>% mutate(cls = map_lgl(rocs, is.null))  %>% dplyr::filter(cls == FALSE)
  #nested <- nested %>% mutate(names = map(rocs, names))

  dumm <- nested %>% dplyr::select(!!sym(config$table$hierarchyKeys()[1]),
                                   !!sym(config$table$hierarchyKeys(TRUE)[1]),
                                   rocs) %>%  unnest()
  dumm <- dumm %>% mutate(comparison = map_chr(rocs, function(x){paste(x$levels, collapse = " ")}))
  dumm <- dumm %>% separate(comparison, into = c("response1" , "response2"), sep=" ")
  dumm <- dumm %>% mutate(auc = map_dbl(rocs, pROC::auc)) %>% arrange(desc(auc))
  return(dumm)
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
#' compute_anova_lm(data, config, hierarchy_level= 2, factor_level=1)
#' compute_anova_lm(data, config, hierarchy_level = 1, factor_level=2)
#' compute_anova_lm(data, config, hierarchy_level= 2, factor_level=2)
compute_anova_lm <- function(data, config, .formula=NULL, hierarchy_level=1, factor_level=1){
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

  pepRes <- data %>% group_by(!!!syms(groupVars)) %>% nest()
  pepRes1 <- pepRes %>% mutate(anova = map( data, aovmodelfit, formula))
  head(pepRes1)
  pepRes2 <- pepRes1 %>% mutate(broomres = map(anova, broom::tidy))
  pepRes3 <- pepRes2 %>%
    dplyr::select(!!!syms(c(groupVars, "broomres"))) %>%
    unnest() %>%
    dplyr::filter(term != "Residuals")

  pVals <- pepRes3 %>% dplyr::select(!!!syms(c(groupVars,"term","p.value")))

  pVals <- pVals %>% mutate(term = glue("{term}.p.value"))
  pVals <- pVals %>% spread("term", "p.value")
  df <- pepRes3 %>% dplyr::select(!!!syms(c(groupVars,"term","df")))
  df <- df %>% mutate(term = glue("{term}.df"))
  df <- df %>% spread(term, df)
  statistic <- pepRes3 %>% dplyr::select(!!!syms(c(groupVars,"term","statistic")))
  statistic <- statistic %>% mutate(term = glue("{term}.statistic"))
  statistic <- statistic %>% spread(term, statistic)
  res <- inner_join(inner_join(pVals, df, by=groupVars), statistic, by=groupVars)
  return(res)
}

# mixed linear models ----

# Creating models from configuration ----

#' get lmer forumula for full model from config
#' @export
model_full_lmer <- function(config, factor_level=2, random= NULL){
  if(factor_level > 2)
  {
    error("can't automatically create model formula")
  }
  formula_str <- paste0(config$table$getWorkIntensity(), " ~ ",
                        "1 + ",
                        paste(config$table$factorKeys()[1:factor_level], collapse=" + "),
                        " + ",
                        paste(config$table$factorKeys()[1:factor_level], collapse=" * "),
                        paste0(" + (1|", config$table$hierarchyKeys(TRUE)[1],")"))
  if(!is.null(random)){
    formula_str <- paste0(formula_str, paste0(" + (1|", random,")"))
  }
  formula <- as.formula( formula_str )

  print(formula)
  res <- function(x, get_formula=FALSE){
    if(get_formula)
    {
      return(formula)
    }

    modelTest <- tryCatch(lmerTest::lmer( formula , data=x ),
                          error=function(e){print(e);return=NULL})
    return(modelTest)
  }
  return(res)
}

#' get mixed model no interactions with peptide random, factor from config
#' @export
model_no_interaction_lmer <- function(config, factor_level=2, random = NULL){
  formula_str <- paste0(config$table$getWorkIntensity(), " ~ ",
                        paste(config$table$factorKeys()[1:factor_level], collapse="+"),
                        paste0(" + (1|", config$table$hierarchyKeys(TRUE)[1],")"))
  if(!is.null(random)){
    formula_str <- paste0(formula_str, paste0(" + (1|", random,")"))
  }
  formula <- as.formula(
    formula_str
  )
  print(formula)
  res <- function(x, get_formula=FALSE){
    if(get_formula)
    {
      return(formula)
    }
    modelTest <- tryCatch(lmerTest::lmer( formula , data=x ),
                          error=function(e){print(e);return=NULL})
    return(modelTest)
  }
  return(res)
}

#' Create custom model
#' @export
#'
make_custom_model <- function(modelstr){
  formula <- as.formula(modelstr)
  res <- function(x, get_formula=FALSE){
    if(get_formula)
    {
      return(formula)
    }
    modelTest <- tryCatch(lmerTest::lmer( formula , data=x ),
                          error=function(e){print(e) ; return=NULL})
    return(modelTest)
  }
  return(res)
}


#' DEPRECATED - get mixed model no interactions with peptide and sample random factor, from config
#' @export
model_no_interaction_and_sample_lmer <- function(config, factor_level = 2){
  warning("DEPRECATED")
  formula <- as.formula(paste0(config$table$getWorkIntensity(), " ~ ",
                               paste(config$table$factorKeys()[1:factor_level], collapse="+"),
                               paste0(" + (1|", config$table$hierarchyKeys(TRUE)[1],")"),
                               paste0(" + (1|", config$table$sampleName,")"),
  ))
  print(formula)
  res <- function(x){
    modelTest <- tryCatch(lmerTest::lmer( formula , data=x ),
                          error=function(e){print(e);return=NULL})
    return(modelTest)
  }
  return(res)
}

#' Likelihood ratio test
#' @export
#'
likelihood_ratio_test <- function(modelNO, model) {
  res <- tryCatch(  anova(modelNO,model), error = function(x) NA)
  if(!is.na(res)){
    res <- broom::tidy(res)[2,"p.value"]
  }
  return(as.numeric(res))
}


# extracting results ----

#' get all comparisons
.contrast_tukey_multcomp <- function(model, factor){
  mcpDef <- multcomp::mcp(Dummy="Tukey")
  names(mcpDef) <- factor
  glt <- multcomp::glht(model, mcpDef)
  sglt <- broom::tidy(summary(glt))
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
  protein_Id <- pepConfig$table$hierarchyKeys()[1]

  for(factor in factors){
    print(factor)
    modelProteinF <- modelProteinF %>%
      mutate(!!paste0("factor_",factor) := purrr::map(!!sym(paste0("lmer_",modelName )), ~.contrast_tukey_multcomp(.,factor=factor)))
  }
  dd <- modelProteinF %>% dplyr::select(protein_Id, starts_with("factor_"))
  contrasts <- dd %>% gather("factor", "contrasts",-!!sym( protein_Id)) %>% unnest() %>% arrange(!!sym(protein_Id)) %>% dplyr::select(-rhs)

  return(contrasts)
}


#' compute all contrasts from non interaction model automatically.
#' @export
#'
workflow_model_contrasts_no_interaction <- function(modelProteinF,
                                                    modelName,
                                                    pepConfig )
{
  result <- list()

  contrasts <- .compute_contrasts_no_interaction(modelProteinF, pepConfig, modelName )
  contrasts <- inner_join(dplyr::select(modelProteinF, protein_Id, isSingular, peptide_Id_n), contrasts )

  result$contrasts <- contrasts
  result$fig <- list()
  result$fig$histogram_coeff_p.values_name <- paste0("Contrasts_Auto_histogram_",modelName,".pdf")
  result$fig$histogram_coeff_p.values <- ggplot(data=contrasts,
                                                aes(x = p.value, group=lhs)) + geom_histogram(bins = 20) + facet_wrap(~lhs)


  result$fig$VolcanoPlot_name <- paste0("Contrasts_Auto_Volcano_",modelName,".pdf")
  result$fig$VolcanoPlot <- quantable::multigroupVolcano(contrasts,
                                                         effect = "estimate",
                                                         type = "p.value",
                                                         condition = "lhs",
                                                         label = "protein_Id",
                                                         xintercept = c(-1, 1),colour = "isSingular")
  result$contrasts <- contrasts
  return(result)
}

#' compute all contrasts from model with interactions based on linfct matrix
#' @export
#'
workflow_model_contrasts_with_interaction <- function(modelProteinF_Int, modelName, linfct ){
  results <- list()

  interaction_model_matrix <- modelProteinF_Int %>%
    mutate(contrasts = map(!!sym(paste0("lmer_",modelName)) , LFQService::my_glht, linfct = linfct, sep=TRUE ))

  interaction_model_matrix %>%
    dplyr::select(protein_Id, contrasts ) %>% unnest() -> modelWithInteractionsContrasts

  modelWithInteractionsContrasts <- inner_join(dplyr::select(modelProteinF_Int, protein_Id, isSingular,peptide_Id_n ),
                                               modelWithInteractionsContrasts)


  results$contrasts <- modelWithInteractionsContrasts
  results$fig <- list()
  results$fig$histogram_coeff_p.values_name <- paste0("Contrasts_histogram_p.values_", modelName ,".pdf")

  results$fig$histogram_coeff_p.values <- ggplot(data=modelWithInteractionsContrasts, aes(x = p.value)) +
    geom_histogram(bins = 20) +
    facet_wrap(~lhs)

  results$fig$VolcanoPlot_name <- paste0("Contrasts_Volcano_",modelName,".pdf")
  results$fig$VolcanoPlot <- quantable::multigroupVolcano(modelWithInteractionsContrasts,
                                                          effect = "estimate",
                                                          type = "p.value",
                                                          condition = "lhs",
                                                          label = "protein_Id",
                                                          xintercept = c(-1, 1), colour = "isSingular")

  return(results)
}


#' helper function to write the result of `workflow_model_contrasts_with_interaction`
#' @export
#' @examples
#' if(FALSE){
#' contrast_interactions <- workflow_model_contrasts_with_interaction(modelProteinF_Int, modelName_Int, XX)
#' modelWithInteractionsContrasts <- inner_join(contrast_interactions$contrasts, likelihood_ratio_test_result )
#' modelWithInteractionsContrasts_Pivot <- pivot_model_contrasts_2_Wide(modelWithInteractionsContrasts)
#' }
write_figures_model_contrasts <- function(contrasts_result,path, fig.width = 10, fig.height = 10){
  pdf(file.path(path,contrasts_result$fig$histogram_coeff_p.values_name), width = fig.width, height = fig.height)
  print(contrasts_result$fig$histogram_coeff_p.values)
  dev.off()

  pdf(file.path(path,contrasts_result$fig$VolcanoPlot_name), width = fig.width, height = fig.height)
  print(contrasts_result$fig$VolcanoPlot)
  dev.off()
}



#' pivot model contrasts matrix to wide format produced by `workflow_model_contrasts_with_interaction` and ...
#' @export
#'
pivot_model_contrasts_2_Wide <- function(modelWithInteractionsContrasts){
  modelWithInteractionsContrasts %>% dplyr::select(protein_Id,likelihood_ratio_test.pValue, lhs, estimate) %>%
    mutate(lhs = glue::glue('estimate.{lhs}')) %>%
    tidyr::spread(lhs, estimate ) -> modelWithInteractionsEstimate

  modelWithInteractionsContrasts %>% dplyr::select(protein_Id,likelihood_ratio_test.pValue, lhs, p.value) %>%
    mutate(lhs = glue::glue('p.value.{lhs}')) %>%
    tidyr::spread(lhs, p.value ) -> modelWithInteractions.p.value

  modelWithInteractionsContrasts_Pivot <- inner_join(modelWithInteractionsEstimate,
                                                     modelWithInteractions.p.value)
  return(modelWithInteractionsContrasts_Pivot)
}

#' get all model coefficients
#' @export
coef_df <-  function(x){
  x <- coef(summary(x));
  x<- data.frame(row.names(x), x);
  return(x)
}

#' run analysis of variance on model and get results
#' @export
anova_df <- function(x){
  x <- anova(x)
  colnames(x) <- make.names(colnames(x))
  x <- data.frame(rownames(x), x)
  return(x)
}

# visualize modelling results ----

#' Plot prdictions
#' @export
#'
plot_lme4_peptide_predictions <- function(m){
  data <- m@frame
  data$prediction <- predict(m)
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  data <- make_interaction_column(data, interactionColumns, sep=":")
  gg <- ggplot(data, aes(x = interaction , y= transformedIntensity)) + geom_point()
  gg <- gg + geom_point(aes(x = interaction, y = prediction), color=2) + facet_wrap(~peptide_Id)
  gg <- gg + theme(axis.text.x=element_text(angle = -90, hjust = 0))
  return(gg)
}


#' plot peptide intensities per interaction with random effects removed
#' @export
#'
plot_lmer4_peptide_noRandom <- function(m,legend.position="none"){
  data <- m@frame
  ran <- lme4::ranef(m)[[1]]
  randeffect <- setdiff(all.vars( terms(formula(m)) ) , all.vars(terms(m)))
  ran <- as.tibble(ran,rownames = randeffect)
  colnames(ran) <- gsub("[()]","",colnames(ran))
  ran <- inner_join(data, ran, by=randeffect)

  ran <- ran %>% mutate(int_randcorrected  = transformedIntensity  - Intercept)
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  ran <- make_interaction_column(ran,interactionColumns, sep=":" )

  meanx <- function(x){mean(x,na.rm=TRUE)}
  gg <- ggplot(ran,aes(x = interaction , y= int_randcorrected, color=peptide_Id)) +
    geom_point(position = position_jitterdodge())
  gg <- gg + stat_summary(fun.y=meanx, colour="black", geom="point",
                          shape=12, size=3,show.legend = FALSE)
  gg <- gg + theme(axis.text.x=element_text(angle = -90, hjust = 0), legend.position =legend.position)
  gg <- gg + geom_boxplot(alpha=0.1)
  return(gg)
}

#' plot intensities per interaction with Two independent random effects removed (1|A) + (1|B)
#' @export
#'
plot_lmer4_peptide_noRandom_TWO <- function(m, legend.position = "none", firstlast = TRUE){

  updateDataWithRandom <- function(data, m, i, randeffect){

    rand_i <- randeffect[i]
    ran <- ranef(m)[[rand_i]]
    name <- paste0(gsub("[()]","",colnames(ran)),"_", rand_i)
    colnames(ran) <- name
    ran <- as.tibble(ran,rownames = rand_i)
    ran <- inner_join(data, ran, by=rand_i)
    ran_res <- ran %>% mutate(int_randcorrected  = transformedIntensity  - !!sym(name))
    ran_res
  }

  data <- m@frame
  if(firstlast){
    i1 <- 1;  i2 <- 2
  } else {
    i1 <- 2;  i2 <- 1
  }
  randeffect <- setdiff(all.vars( terms(formula(m)) ) , all.vars(terms(m)))
  ran_res <- updateDataWithRandom(data, m, i1, randeffect)
  ran_res <- updateDataWithRandom(ran_res, m, i2, randeffect)

  #randeffect <- setdiff(all.vars( terms(formula(m)) ) , all.vars(terms(m)))
  #rand_i <- randeffect[i]
  #ran <- ranef(m)[[rand_i]]
  #name <- paste0(gsub("[()]","",colnames(ran)),"_", rand_i)
  #colnames(ran) <- name
  #ran <- as.tibble(ran,rownames = rand_i)
  #ran_res <- inner_join(ran_res, ran, by=rand_i)
  #ran_res <- ran_res %>% mutate(int_randcorrected  = int_randcorrected  - !!sym(name))
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  ran_res <- make_interaction_column(ran_res,interactionColumns, sep=":" )

  meanx <- function(x){mean(x,na.rm=TRUE)}
  gg <- ggplot(ran_res,aes(x = interaction , y= int_randcorrected, color=!!sym(randeffect[i1]))) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2))
  gg <- gg + stat_summary(fun.y=meanx, colour="black", geom="point",
                          shape=12, size=3,show.legend = FALSE)
  gg <- gg + theme(axis.text.x=element_text(angle = -90, hjust = 0),legend.position=legend.position)
  gg <- gg + geom_boxplot(alpha=0.1)
  gg
  return(gg)
}


#' Add predicted values for each interaction
#' @export
#'
plot_predicted_interactions <- function(gg, m){
  cm <- .lmer4_coeff_matrix(m)
  xstart_end <- data.frame(xstart = rownames(cm$mm), xend = rownames(cm$mm))
  ystart_end <- data.frame(xend = rownames(cm$mm), ystart =rep(0, nrow(cm$mm)),
                           yend = cm$mm %*% cm$coeffs)
  segments <- inner_join(xstart_end, ystart_end, by="xend")
  gg <- gg + geom_segment(aes(x = xstart, y = ystart , xend = xend, yend =yend), data=segments, color = "blue", arrow=arrow())
  return(gg)
}

#' Make model plot with title - protein Name.
#' @export
#'
plot_model_and_data <- function(m, proteinID, legend.position = "none"){
  gg <- plot_lmer4_peptide_noRandom(m,legend.position=legend.position)
  gg <- plot_predicted_interactions(gg, m)
  gg <- gg + ggtitle(proteinID)
  gg
}


#' Plotting two independent random effects
#' @export
#'
plot_model_and_data_TWO <- function(m, proteinID, legend.position = "none" , firstlast= TRUE){
  gg <- plot_lmer4_peptide_noRandom_TWO(m, legend.position = legend.position, firstlast = firstlast)
  gg <- plot_predicted_interactions(gg, m)
  gg <- gg + ggtitle(proteinID)
  gg
}


# generate linear functions -----

#' get matrix of indicator coefficients for each interaction
#'
.lmer4_coeff_matrix <- function(m){
  data <- m@frame
  interactionColumns <- intersect(attributes(terms(m))$term.labels,colnames(data))
  data <- make_interaction_column(data, interactionColumns, sep=":")

  coeffs <- coefficients(summary(m))[,'Estimate']

  inter <- unique(data$interaction)
  mm <- matrix(0, nrow=length(inter), ncol=length(coeffs))
  rownames(mm) <- inter
  colnames(mm) <- names(coeffs)
  mm[,1]<-1
  coefi <- coeffs[-1]
  for(i in 1:length(coefi)){
    positionIDX <- grep(names(coefi)[i], inter)
    mm[positionIDX,i+1] <- 1
  }
  return(list(mm = mm, coeffs = coeffs))
}


#' coeff_weights_factor_levels
.coeff_weights_factor_levels <- function(mm){
  getCoeffs <- function(factor_level, mm){
    idx <- grep(factor_level, rownames(mm))
    x <- as.list(apply(mm[idx,, drop=FALSE],2,mean) )
    x <- as.tibble(x)
    add_column(x, "factor_level" = factor_level,.before=1)
  }
  factor_levels <- unique(unlist(str_split(rownames(mm), ":")))
  xx <- map_df(factor_levels, getCoeffs, mm)
  return(xx)
}

#' get linfct from model
#' @export
#' @examples
#' if(FALSE){
#' m <- modelProteinF_Int$lmer_visittime_group_interactions[[1]]
#' linfct <- lmer4_linfct_from_model(m)
#' linfct$linfct_factors
#' linfct$linfct_interactions
#' }
#'
lmer4_linfct_from_model <- function(m){
  cm <- .lmer4_coeff_matrix(m)
  cm_mm <- cm$mm[order(rownames(cm$mm)),]

  dd <- .coeff_weights_factor_levels(cm_mm)
  dd_m <- dd %>% dplyr::select(-factor_level) %>% data.matrix()
  rownames(dd_m) <- dd$factor_level
  dd_m <- dd_m[order(rownames(dd_m)),]

  return(list(linfct_factors = dd_m , linfct_interactions = cm_mm))
}



# Computing contrasts helpers -----

#' apply glht method to linfct
#' @export
my_glht <- function(model , linfct , sep=FALSE){
  if(sep){
    res <- list()
    for(i in 1:nrow(linfct)){
      x <- glht(model, linfct=linfct[i,,drop=FALSE])
      RHS <- broom::tidy(confint(x)) %>% dplyr::select(-estimate)
      x <- inner_join(broom::tidy(summary(x)),RHS,by = c("lhs", "rhs")) %>% dplyr::select(-rhs)
      res[[i]] <- x
    }
    res <- bind_rows(res)
    return(res)
  }else{
    x <- glht(model, linfct = linfct)
    RHS <- broom::tidy(confint(x)) %>% dplyr::select(-estimate)
    res <- inner_join(broom::tidy(summary(x)),RHS,by = c("lhs", "rhs")) %>% dplyr::select(-rhs)
    res
  }
}



#' Workflow function to apply modelFunction to peroteinList nestProtein
#'
#' apply modelling, extracts coefficients,
#' funs anova, filters results, generates histogram of p-values, pairsplot
#' @export
#' @examples
#' allData_PhonixDS_1097969 <- readRDS("data/allData_PhonixDS_1097969.Rds")
#' D <- allData_PhonixDS_1097969
#' names(D)
#' usethis::use_data(allData_PhonixDS_1097969)

workflow_lme4_model_analyse <- function(nestProtein, modelFunction, modelName, prot_stats)
{
  lmermodel <- paste0("lmer_", modelName)
  exists_lmer <- paste0("exists_lmer_", modelName)
  Coeffs_model <- paste0("Coeffs_", modelName)
  Anova_model <- paste0("Anova_", modelName)

  nestProtein %>% mutate(!!lmermodel := purrr::map(data, modelFunction)) -> modelProtein

  modelProtein <- modelProtein %>% mutate(!!exists_lmer := map_lgl(!!sym(lmermodel), function(x){!is.null(x)}))
  modelProteinF <- modelProtein %>% dplyr::filter( !!sym(exists_lmer) == TRUE)
  no_ModelProtein <- modelProtein %>% dplyr::filter(!!sym(exists_lmer) == FALSE)

  modelProteinF <- modelProteinF %>% mutate(isSingular = map_lgl(!!sym(lmermodel), lme4::isSingular ))
  modelProteinF <- modelProteinF %>% mutate(nrcoef = map_int(!!sym(lmermodel), function(x){ncol(coef(x)[[1]])} ))

  modelProteinF <- modelProteinF %>% mutate(!!Coeffs_model := purrr::map(!!sym(lmermodel),  coef_df ))
  modelProteinF <- modelProteinF %>% mutate(!!Anova_model := purrr::map(!!sym(lmermodel),  anova_df ))
  modelProteinF <- inner_join(modelProteinF, prot_stats)

  Model_Coeff <- modelProteinF %>% dplyr::select(protein_Id, !!sym(Coeffs_model), isSingular, nrcoef) %>% unnest()
  Model_Anova <- modelProteinF %>% dplyr::select(protein_Id, !!sym(Anova_model), isSingular, nrcoef) %>% unnest()

  reslist <- list()
  reslist$models <- modelProteinF
  reslist$no_models <- no_ModelProtein

  reslist$fig <- list()
  reslist$table <- list()

  reslist$fig$histogram_coeff_p.values <- ggplot(data=Model_Coeff, aes(x = Pr...t.., group=row.names.x.)) +
    geom_histogram(bins = 20) +
    facet_wrap(~row.names.x.)

  reslist$fig$VolcanoPlot <- Model_Coeff %>% dplyr::filter(row.names.x. != "(Intercept)") %>%  quantable::multigroupVolcano(
    effect = "Estimate",
    type = "Pr...t..",
    condition = "row.names.x.",
    label = "protein_Id", xintercept = c(-1, 1),colour = "isSingular")

  Model_Coeff %>% dplyr::select(protein_Id , row.names.x. ,  Estimate ) %>% tidyr::spread(row.names.x.,Estimate ) -> forPairs
  reslist$fig$Pairsplot_Coef <-  GGally::ggpairs(forPairs, columns=2:ncol(forPairs))


  reslist$fig$histogram_anova_p.values <- ggplot(data=Model_Anova, aes(x = Pr..F., group=rownames.x.)) +
    geom_histogram(bins = 20) +
    facet_wrap(~rownames.x.)

  reslist$table$Model_Coeff <- inner_join(prot_stats, Model_Coeff)
  reslist$table$Model_Anova <-inner_join( prot_stats , Model_Anova )

  return(reslist)
}

#' Writes figures generated by `workflow_lme4_model_analyse`
#' @export
write_figures_lme4_model_analyse <- function(modelling_result, modelName, path, fig.width = 10 , fig.height = 10){
  pdf(file.path(path,paste0("Coef_Histogram_",modelName,".pdf")), width = fig.width, height = fig.height )
  print(modelling_result$fig$histogram_coeff_p.values)
  dev.off()

  pdf(file.path(path,paste0("Coef_VolcanoPlot_",modelName,".pdf")), width = fig.width , height = fig.height)
  print(modelling_result$fig$VolcanoPlot)
  dev.off()

  pdf(file.path(path, paste0("Coef_Pairsplot_",modelName,".pdf")), width = fig.width , height = fig.height)
  print(modelling_result$fig$Pairsplot_Coef)
  dev.off()

  pdf(file.path(path,paste0("Anova_p.values_", modelName, ".pdf")), width = fig.width , height = fig.height)
  histogram_anova_p.values <-
    print(modelling_result$fig$histogram_anova_p.values)
  dev.off()
}


#' p2621 workflow interaction
#' @export
#' @examples
#'
#' D <- LFQService::allData_PhonixDS_1097969
#' modelName  <- "f_Condition_r_peptide"
#' formula_randomPeptide <- make_custom_model("transformedIntensity  ~ Condition + (1 | peptide_Id)")
#' res_cond_r_pep <- workflow_no_interaction_modelling(D,formula_randomPeptide, modelName)
workflow_interaction_modelling <- function(results, modelFunction, modelName, writeCoefResults=FALSE){

  pepIntensity <- results$pepIntensityNormalized
  pepConfig <- results$config_pepIntensityNormalized

  pepIntensity %>%
    group_by(!!sym(pepConfig$table$hierarchyKeys()[1])) %>%
    nest() -> nestProtein

  prot_stats <- summarizeHierarchy(pepIntensity, pepConfig)

  interactionResults <- workflow_lme4_model_analyse(nestProtein, modelFunction, modelName, prot_stats = prot_stats)

  if(writeCoefResults){
    write_figures_lme4_model_analyse(interactionResults, modelName, results$path)
    readr::write_csv(interactionResults$table$Model_Coeff,
                     path = file.path( results$path, paste0("Coef_",modelName, ".txt")))
    readr::write_csv(interactionResults$table$Model_Anova,
                     path = file.path( results$path , paste0("ANOVA_",modelName,".txt" ) ))
  }

  modelProteinF_Int <- interactionResults$models
  modelProteinF_Int <- modelProteinF_Int %>%
    dplyr::filter(nrcoef==as.numeric(tail(names(table(modelProteinF_Int$nrcoef)),n=1)))

  return(modelProteinF_Int)
}

#' p2621 workflow no interaction - adds contrast computation to `workflow_interaction_modelling` function.
#' @export
workflow_no_interaction_modelling <- function(results, modelFunction, modelName,writeCoefResults=FALSE){
  modelProteinF <- workflow_interaction_modelling(results, modelFunction, modelName, writeCoefResults=FALSE)

  contrasts <- workflow_model_contrasts_no_interaction(modelProteinF,
                                                       modelName ,
                                                       pepConfig )
  write_figures_model_contrasts(contrasts, results$path)

  return(list(TwoFactorModelFactor2 = contrasts$contrasts, modelProteinF = modelProteinF, modelName = modelName))
}


#' p2621 workflow likelihood ratio test
#' @export
workflow_likelihood_ratio_test <- function(modelProteinF,modelName, modelProteinF_Int,modelName_Int ){
  # Model Comparison
  reg <- inner_join(dplyr::select(modelProteinF, protein_Id, starts_with("lmer_")),
                    dplyr::select(modelProteinF_Int, protein_Id, starts_with("lmer_")) , by="protein_Id")

  reg <- reg %>% mutate(modelComparisonLikelihoodRatioTest = map2(!!sym(paste0("lmer_", modelName)),
                                                                  !!sym(paste0("lmer_", modelName_Int)),
                                                                  likelihood_ratio_test ))
  likelihood_ratio_test_result <- reg %>%
    dplyr::select(protein_Id, modelComparisonLikelihoodRatioTest) %>% unnest()
  likelihood_ratio_test_result <- likelihood_ratio_test_result %>%
    dplyr::rename(likelihood_ratio_test.pValue = modelComparisonLikelihoodRatioTest)
  return(likelihood_ratio_test_result)
}

#' p2621 compute group averages
#' @export
workflow_group_averages <- function(models,
                                    modelName,
                                    path,
                                    lin_int,
                                    likelihood_ratio_test_result = NULL ){
  #computeGroupAverages
  interaction_model_matrix <- models %>%
    mutate(groupAverages = map(!!sym(paste0("lmer_", modelName)) , LFQService::my_glht, linfct = lin_int, sep=TRUE ))

  interaction_model_matrix %>%
    dplyr::select(protein_Id, groupAverages ) %>% unnest() -> groupAverages

  if(!is.null(likelihood_ratio_test_result)){
    groupAverages <- inner_join(groupAverages,  likelihood_ratio_test_result)
  }

  groupAveragesW <- groupAverages %>% dplyr::select(-p.value, - statistic)

  #  separate(protein_Id, c("type", "Uniprot", "ProteinName"), sep="\\|", remove = FALSE)
  write_csv(groupAveragesW, path=file.path(path, paste0("Group_Averages_", modelName, ".csv")))
  return(groupAverages)
}


#' p2621 workflow linfunct contrasts
#' @export
workflow_linfunct_contrasts <- function(models, modelName, likelihood_ratio_test_result ,  linfct, path,
                                        fig.width=10, fig.height=10)
{
  contrast_interactions <- workflow_model_contrasts_with_interaction(models, modelName, linfct)
  write_figures_model_contrasts(contrast_interactions, path,
                                fig.width = fig.width , fig.height = fig.height)

  modelWithInteractionsContrasts <- inner_join(contrast_interactions$contrasts,
                                               likelihood_ratio_test_result )
  write_csv(modelWithInteractionsContrasts,
            path=file.path(path, paste0("Contrasts_SignificanceValues_", modelName, ".csv")))

  modelWithInteractionsContrasts_Pivot <- pivot_model_contrasts_2_Wide(modelWithInteractionsContrasts)
  write_csv(modelWithInteractionsContrasts_Pivot,
            path=file.path(path, paste0("Contrasts_SignificanceValues_", modelName, "_PIVOT.csv")))

  return(list(contrast_interactions = contrast_interactions,
              modelWithInteractionsContrasts = modelWithInteractionsContrasts,
              modelWithInteractionsContrasts_Pivot= modelWithInteractionsContrasts_Pivot
              ))
}

