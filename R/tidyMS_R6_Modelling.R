# function for modelling go here.
#' rocs helper function
rocs <- function(data ,response, predictor){
  responseX <- data %>% dplyr::pull(!!sym(response))
  predictorX <- data %>% dplyr::pull(!!sym(predictor))
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
#'
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
#'
compute_roc <- function(data, config){
  nested <- data %>% group_by(!!sym(config$table$hierarchyKeys()[1]) ,
                              !!sym(config$table$hierarchyKeys(TRUE)[1])) %>% nest()
  nested <- nested %>% dplyr::mutate(rocs = map(data ,
                                         rocs, response = config$table$factorKeys()[1],
                                         predictor= config$table$getWorkIntensity() ))

  nested <- nested %>% dplyr::mutate(cls = map_lgl(rocs, is.null))  %>%
    dplyr::filter(cls == FALSE)
  #nested <- nested %>% mutate(names = map(rocs, names))

  dumm <- nested %>% dplyr::select(!!sym(config$table$hierarchyKeys()[1]),
                                   !!sym(config$table$hierarchyKeys(TRUE)[1]),
                                   rocs) %>%  tidyr::unnest()
  dumm <- dumm %>% dplyr::mutate(comparison = map_chr(rocs, function(x){paste(x$levels, collapse = " ")}))
  dumm <- dumm %>% separate(comparison, into = c("response1" , "response2"), sep=" ")
  dumm <- dumm %>% dplyr::mutate(auc = map_dbl(rocs, pROC::auc)) %>%
    arrange(desc(auc))
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

  pepRes <- data %>% dplyr::group_by(!!!syms(groupVars)) %>% nest()
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

#' Create custom lmer model
#' @export
#'
make_custom_model_lmer <- function( modelstr ) {
  formula <- as.formula(modelstr)
  res <- function(x, get_formula=FALSE){
    if(get_formula)
    {
      return(formula)
    }
    modelTest <- tryCatch(lmerTest::lmer( formula , data=x ),
                          error = function(e){print(e) ; return=NULL})
    return(modelTest)
  }
  return(res)
}

#' Create custom ml model
#' @export
#'
make_custom_model_lm <- function( modelstr ) {
  formula <- as.formula(modelstr)
  res <- function(x, get_formula=FALSE){
    if(get_formula)
    {
      return(formula)
    }
    modelTest <- tryCatch(lm( formula , data=x ),
                          error = function(e){print(e) ; return=NULL})
    return(modelTest)
  }
  return(res)
}


.likelihood_ratio_test <- function(modelNO, model) {
  res <- tryCatch(  anova(modelNO,model), error = function(x) NULL)
  if(!is.null(res)){
    res <- suppressWarnings(broom::tidy(res))[2,"p.value"]
    return(as.numeric(res))

  }
  else{
    return(NA)
  }
}


# computing contrast ----

#' compute group averages
#'
#' used in p2621, p2109
#'
#' @export
workflow_contrasts_linfct <- function(models,
                                      modelName,
                                      linfct,
                                      subject_Id = "protein_Id"){
  #computeGroupAverages
  modelcol <-paste0("lmer_", modelName)
  contrastfun <- NULL
  whatclass <- class(models[[modelcol]][[1]])


  if(whatclass=="lm"){
    contrastfun = LFQService::my_glht
  }else if(whatclass == "lmerModLmerTest"){
    contrastfun = LFQService::my_contest
  }

  interaction_model_matrix <- models %>%
    dplyr::mutate(contrast = map(!!sym(modelcol) , contrastfun , linfct = linfct, sep=TRUE ))
  mclass <- function(x){
    class(x)[1]
  }
  interaction_model_matrix %>%
    dplyr::mutate(classC = map_chr(contrast,mclass)) %>%
    filter(classC != "logical") -> interaction_model_matrix

  interaction_model_matrix %>%
    dplyr::select_at( c(subject_Id, "contrast") ) %>% tidyr::unnest() -> contrasts

  isSing <- models %>% dplyr::select_at(c(subject_Id, "isSingular")) %>% distinct()
  contrasts <- inner_join(contrasts, isSing, by=subject_Id)
  return(contrasts)
}

#' write results of `workflow_contrasts_linfct`
#' @export
#'
workflow_contrasts_linfct_write <- function(results, modelName, path, prefix = "Contrasts", subject_Id = "protein_Id"){
  fileLong <-file.path(path,paste0(prefix,"_",modelName,".csv"))
  readr::write_csv(results, path = fileLong)
  fileWide <-file.path(path,paste0(prefix,"_",modelName,"_PIVOT.csv"))
  resultswide <- pivot_model_contrasts_2_Wide(results, subject_Id = subject_Id)
  readr::write_csv(resultswide, path = fileWide)

}



#' visualize output of `workflow_contrasts_linfct``
#' @export
#'
workflow_contrasts_linfct_vis <- function(contrasts, modelName, prefix = "Contrasts"){
  fig <- list()
  fig$histogram_coeff_p.values_name <- paste0(prefix,"_Histogram_p.values_", modelName ,".pdf")

  fig$histogram_coeff_p.values <- ggplot(data=contrasts, aes(x = p.value)) +
    geom_histogram(bins = 20) +
    facet_wrap(~lhs)

  fig$VolcanoPlot_name <- paste0(prefix,"_Volcano_",modelName,".pdf")
  fig$VolcanoPlot <- quantable::multigroupVolcano(contrasts,
                                                  effect = "estimate",
                                                  type = "p.value",
                                                  condition = "lhs",
                                                  label = "protein_Id",
                                                  xintercept = c(-1, 1), colour = "isSingular")

  return(fig)
}


#' helper function to write the result of `workflow_contrasts_linfct_vis`
#'
#' used in p2901
#'
#' @export
workflow_contrasts_linfct_vis_write <- function(contrasts_result, path, fig.width = 10, fig.height = 10){
  pdf(file.path(path,contrasts_result$histogram_coeff_p.values_name), width = fig.width, height = fig.height)
  print(contrasts_result$histogram_coeff_p.values)
  dev.off()

  pdf(file.path(path,contrasts_result$VolcanoPlot_name), width = fig.width, height = fig.height)
  print(contrasts_result$VolcanoPlot)
  dev.off()
}

#' pivot model contrasts matrix to wide format produced by `workflow_contrasts_linfct` and ...
#' @export
#'
pivot_model_contrasts_2_Wide <- function(modelWithInteractionsContrasts, subject_Id = "protein_Id"){
  modelWithInteractionsContrasts %>%
    dplyr::select_at(c(subject_Id, "lhs", "estimate")) %>%
    dplyr::mutate(lhs = glue::glue('estimate.{lhs}')) %>%
    tidyr::spread(lhs, estimate ) -> modelWithInteractionsEstimate

  modelWithInteractionsContrasts %>%
    dplyr::select_at(c(subject_Id, "lhs", "p.value")) %>%
    dplyr::mutate(lhs = glue::glue('p.value.{lhs}')) %>%
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

# visualize lmer modelling results ----

#' Plot prdictions
#' @export
#' @examples
#' m <- LFQService::interactionModel_p1807
#' plot_lmer_peptide_predictions(m)
plot_lmer_peptide_predictions <- function(m){
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
#' @examples
#' library(tidyverse)
#' m <-LFQService::basicModel_p1807
#' plot_lmer_peptide_noRandom(m)
#'
#' m <- LFQService::interactionModel_p1807
#' plot_lmer_peptide_noRandom(m)
plot_lmer_peptide_noRandom <- function(m,legend.position="none"){
  data <- m@frame
  ran <- lme4::ranef(m)[[1]]
  randeffect <- setdiff(all.vars( terms(formula(m)) ) , all.vars(terms(m)))
  ran <- tibble::as_tibble(ran,rownames = randeffect)
  colnames(ran) <- gsub("[()]","",colnames(ran))
  ran <- inner_join(data, ran, by=randeffect)

  ran <- ran %>% dplyr::mutate(int_randcorrected  = transformedIntensity  - Intercept)
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

#' plot intensities per interaction with Two
#'  independent random effects removed (1|A) + (1|B)
#' @export
#' @examples
#' m <- LFQService::interactionModel_p1807
#' #plot_lmer_peptide_noRandom_TWO(m)
plot_lmer_peptide_noRandom_TWO <- function(m, legend.position = "none", firstlast = TRUE){

  updateDataWithRandom <- function(data, m, i, randeffect){
    rand_i <- randeffect[i]
    ran <- lme4::ranef(m)[[rand_i]]
    name <- paste0(gsub("[()]","",colnames(ran)),"_", rand_i)
    colnames(ran) <- name
    ran <- tibble::as_tibble(ran,rownames = rand_i)
    ran <- inner_join(data, ran, by=rand_i)
    ran_res <- ran %>% dplyr::mutate(int_randcorrected  = transformedIntensity  - !!sym(name))
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
  #ran <- as_tibble(ran,rownames = rand_i)
  #ran_res <- inner_join(ran_res, ran, by=rand_i)
  #ran_res <- ran_res %>% dplyr::mutate(int_randcorrected  = int_randcorrected  - !!sym(name))
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
#' @examples
#' #m <- LFQService::interactionModel_p1807
#' #plot_lmer_predicted_interactions(m)
plot_lmer_predicted_interactions <- function(gg, m){
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
#' @examples
#' m <- LFQService::interactionModel_p1807
#' plot_lmer_model_and_data(m,"dumm")
plot_lmer_model_and_data <- function(m, proteinID, legend.position = "none"){
  gg <- plot_lmer_peptide_noRandom(m,legend.position=legend.position)
  gg <- plot_lmer_predicted_interactions(gg, m)
  gg <- gg + ggtitle(proteinID)
  gg
}


#' Plotting two independent random effects
#' @export
#'
plot_lmer_model_and_data_TWO <- function(m, proteinID, legend.position = "none" , firstlast= TRUE){
  gg <- plot_lmer_peptide_noRandom_TWO(m, legend.position = legend.position, firstlast = firstlast)
  gg <- plot_lmer_predicted_interactions(gg, m)
  gg <- gg + ggtitle(proteinID)
  gg
}


# Generate linear functions -----

#' get matrix of indicator coefficients for each interaction
#'
.lmer4_coeff_matrix <- function(m){
  data <- NULL
  if(class(m)  == "lm"){
    data <- m$model
  }else{
    # for "lmerModLmerTest"
    data <- m@frame
  }
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
    x <- tibble::as_tibble(x)
    tibble::add_column(x, "factor_level" = factor_level,.before=1)
  }
  factor_levels <- unique(unlist(stringr::str_split(rownames(mm), ":")))
  xx <- purrr::map_df(factor_levels, getCoeffs, mm)
  return(xx)
}

#' get linfct from model
#' @export
#' @examples
#' #if(FALSE){
#'
#' m <- LFQService::basicModel_p1807
#' m
#' linfct <- linfct_from_model(m)
#'
#' linfct$linfct_factors
#' linfct$linfct_interactions
#'
#' m <- LFQService::interactionModel_p1807
#' linfct <- linfct_from_model(m)
#' all.equal(linfct$linfct_factors["CelltypeCMP/MEP",] , apply(linfct$linfct_interactions[grep("CelltypeCMP/MEP", rownames(linfct$linfct_interactions)),],2, mean))
#' linfct$linfct_interactions
#' #}
#'
#' m <- lm(Petal.Width ~ Species, data=iris)
#' linfct_from_model(m)
linfct_from_model <- function(m){

  cm <- .lmer4_coeff_matrix(m)
  cm_mm <- cm$mm[order(rownames(cm$mm)),]

  dd <- .coeff_weights_factor_levels(cm_mm)
  dd_m <- dd %>% dplyr::select(-factor_level) %>% data.matrix()
  rownames(dd_m) <- dd$factor_level
  dd_m <- dd_m[order(rownames(dd_m)),]

  return(list(linfct_factors = dd_m , linfct_interactions = cm_mm))
}

#' create all possible contrasts
#' @export
#' @examples
#' m <- LFQService::basicModel_p1807
#' m
#' linfct <- linfct_from_model(m)
#'
#' xl <- linfct_all_possible_contrasts(linfct$linfct_factors)
#' xx <- linfct_all_possible_contrasts(linfct$linfct_interactions)
#'
linfct_all_possible_contrasts <- function( lin_int ){
  combs <- combn(nrow(lin_int),2)
  names <- rownames(lin_int)
  newnames <- rep("", ncol(combs))
  new_lin_fct <- matrix(NA,  nrow= ncol(combs), ncol =ncol(lin_int))
  for(i in 1:ncol(combs)){
    newnames[i] <- paste(names[combs[,i]], collapse=" - ")
    new_lin_fct[i,] <- lin_int[combs[1,i],] - lin_int[combs[2,i],]
  }
  rownames(new_lin_fct) <- newnames
  colnames(new_lin_fct) <- colnames(lin_int)
  return(new_lin_fct)
}
#' create contrasts between factor levels
#' @export
#' @examples
#' library(LFQService)
#' m <- LFQService::basicModel_p1807
#' xl <- linfct_factors_contrasts(m)
#' xl
#' m <- lm(Petal.Width ~ Species, data=iris)
#' linfct_factors_contrasts(m)
linfct_factors_contrasts <- function(m){
  ffac <- attributes(terms(m))$term.labels
  linfct_factors <- linfct_from_model(m)$linfct_factors
  factorLevels <- rownames(linfct_factors)
  res <- vector(length(ffac), mode = "list")
  for(i in 1:length(ffac)){
    fac <- ffac[i]
    idx <- grep(fac, factorLevels)
    res[[i]]<-linfct_all_possible_contrasts(linfct_factors[idx,])
  }
  res <- do.call(rbind, res)
  return(res)
}

# Computing contrasts helpers -----




#' apply multcomp::glht method to linfct
#' @export
#' @examples
#'
#' mb <- LFQService::basicModel_p1807
#' linfct <- linfct_from_model(mb)
#' names(linfct)
#' my_glht(mb, linfct$linfct_factors)
#' my_glht(mb, linfct$linfct_interactions)
#' mi <-  LFQService::interactionModel_p1807
#' linfct_int <- linfct_from_model(mb)
#' names(linfct_int)
#' my_glht(mi, linfct_int$linfct_factors)
#' my_glht(mi, linfct_int$linfct_interactions)
#'
my_glht <- function(model , linfct , sep=FALSE ) {
  if(!class(model) == "lm") # fixes issue of mutlcomp not working on factors of class character
  {
    warning("USE ONLY WITH LM models ", class(model))
    if(length(lme4::fixef(model)) != ncol(linfct) ){
      return(NA) # catch rank defficient
    }
  }else{
    model$model <- as.data.frame(unclass(model$model))
  }
  if(sep){
    res <- list()
    for(i in 1:nrow(linfct)){
      x <- multcomp::glht(model, linfct=linfct[i,,drop=FALSE])
      RHS <- broom::tidy(confint(x)) %>% dplyr::select(-estimate)
      x <- inner_join(broom::tidy(summary(x)),RHS,by = c("lhs", "rhs")) %>% dplyr::select(-rhs)
      res[[i]] <- x
    }
    res <- bind_rows(res)
    return(res)
  }else{
    x <- multcomp::glht(model, linfct = linfct)
    RHS <- broom::tidy(confint(x)) %>% dplyr::select(-estimate)
    res <- inner_join(broom::tidy(summary(x)), RHS, by = c("lhs", "rhs")) %>%
      dplyr::select(-rhs)
    res
  }
}

#' applies contrast computation using lmerTest::contest function
#' @export
#' @examples
#' mb <- LFQService::basicModel_p1807
#' linfct <- linfct_from_model(mb)
#' names(linfct)
#' lmerTest::contest(mb, linfct$linfct_interactions, joint = FALSE, confint = TRUE)
#' my_contest(mb, linfct$linfct_factors)
#' my_contest(mb, linfct$linfct_interactions)
#' my_glht(mb, linfct$linfct_factors)
#' my_glht(mb, linfct$linfct_interactions)
my_contest <- function(model, linfct , sep=TRUE){
  if(length(lme4::fixef(model)) != ncol(linfct) ){
    return(NA) # catch rank defficient
  }
  res <- lmerTest::contest(model, linfct, joint = FALSE, confint = TRUE)
  res <- as_tibble(res, rownames="lhs")
  res <- res %>% rename(estimate = Estimate, std.error = "Std. Error",
                        statistic="t value",
                        p.value = "Pr(>|t|)",
                        conf.low = "lower",
                        conf.high = "upper")
  return(res)
}

# Fit the models to data ----

#' check if lm model is singular
#' @export
#'
isSingular_lm <- function(m){
  anyNA <- any(is.na(coef(m)))
  if(anyNA){
    return(TRUE)
  }else{
    return(FALSE)
  }

}


#' analyses lmer4 and lm models created using help function `make_custom_model_lm` or `make_custom_model_lmer`
#'
#' used in project p2901
#'
#' @export
#' @examples
#' D <- LFQService::resultsV12954
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <- make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)")
#' pepIntensity <- D$pepIntensityNormalized
#' config <- D$config_pepIntensityNormalized
#' config$table$hkeysLevel()
#' modellingResult <- workflow_model_analyse( pepIntensity,
#'  config,
#'  formula_randomPeptide,
#'  modelName,
#'  isSingular = lme4::isSingular)
#'
workflow_model_analyse <- function(pepIntensity,config,
                                   modelFunction, modelName, isSingular = isSingular_lm)
{

  hierarchyKey <- config$table$hkeysLevel()
  pepIntensity %>%
    group_by(!!!syms(hierarchyKey)) %>%
    nest() -> nestProtein
  prot_stats <- summarizeHierarchy(pepIntensity, config)

  lmermodel <- paste0("lmer_", modelName)
  exists_lmer <- paste0("exists_lmer_", modelName)
  Coeffs_model <- paste0("Coeffs_", modelName)
  Anova_model <- paste0("Anova_", modelName)

  nestProtein %>% dplyr::mutate(!!lmermodel := purrr::map(data, modelFunction)) -> modelProtein

  modelProtein <- modelProtein %>% dplyr::mutate(!!exists_lmer := map_lgl(!!sym(lmermodel), function(x){!is.null(x)}))
  modelProteinF <- modelProtein %>% dplyr::filter( !!sym(exists_lmer) == TRUE)
  no_ModelProtein <- modelProtein %>% dplyr::filter(!!sym(exists_lmer) == FALSE)

  modelProteinF <- modelProteinF %>% dplyr::mutate(!!"isSingular" := map_lgl(!!sym(lmermodel), isSingular ))
  nrcoeff <- function(x){
    cc <- coef(x)
    if(class(cc) == "numeric"){
      return(length(cc))
    }else{
      return(ncol(cc[[1]]))
    }
  }

  modelProteinF <- modelProteinF %>% dplyr::mutate(nrcoef = map_int(!!sym(lmermodel), nrcoeff))
  modelProteinF <- modelProteinF %>% dplyr::mutate(!!Coeffs_model := purrr::map( !!sym(lmermodel),  coef_df ))
  modelProteinF <- modelProteinF %>% dplyr::mutate(!!Anova_model := purrr::map( !!sym(lmermodel),  anova_df ))
  modelProteinF <- inner_join(modelProteinF, prot_stats)

  Model_Coeff <- modelProteinF %>% dplyr::select(!!sym(hierarchyKey), !!sym(Coeffs_model), isSingular, nrcoef) %>%
    tidyr::unnest()
  Model_Coeff <- inner_join(prot_stats, Model_Coeff)
  Model_Anova <- modelProteinF %>% dplyr::select(!!sym(hierarchyKey), !!sym(Anova_model), isSingular, nrcoef) %>%
    tidyr::unnest()
  Model_Anova <-inner_join( prot_stats , Model_Anova )

  return(list(modelProteinF = modelProteinF,
              no_ModelProtein = no_ModelProtein,
              modelName = modelName,
              Model_Coeff=Model_Coeff,
              Model_Anova= Model_Anova))
}

#' writes results of `workflow_model_analyse`, anova table and all the coefficients with parameters.
#' @export
workflow_model_analyse_write <- function(modellingResult, modelName, path){

  message("writing tables into :", path)
  readr::write_csv(modellingResult$Model_Coeff,
                   path = file.path( path, paste0("Coef_",modelName, ".txt")))
  readr::write_csv(modellingResult$Model_Anova,
                   path = file.path( path , paste0("ANOVA_",modelName,".txt" ) ))
}

#' visualize workflow model analyse results
#'
#' used in p2901
#'
#' @export
#' @examples
#'
#' D <- LFQService::resultsV12954
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <- make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)")
#' modellingResult <-  workflow_model_analyse(D$pepIntensityNormalized,D$config_pepIntensityNormalized, formula_randomPeptide, modelName)
#' res <- workflow_model_analyse_vis(modellingResult,D$config_pepIntensityNormalized, modelName)
#' #workflow_model_analyse_vis_write(res,modelName,path = ".")
#'
workflow_model_analyse_vis <- function(modellingResult, config, modelName) {

  lmermodel <- paste0("lmer_", modelName)
  exists_lmer <- paste0("exists_lmer_", modelName)
  Coeffs_model <- paste0("Coeffs_", modelName)
  Anova_model <- paste0("Anova_", modelName)
  Model_Coeff <- modellingResult$Model_Coeff
  Model_Anova <- modellingResult$Model_Anova

  hierarchyKey <- config$table$hierarchyKeys()[1]
  modelProteinF <- modellingResult$modelProteinF

  fig <- list()

  fig$histogram_coeff_p.values <- ggplot(data = Model_Coeff, aes(x = Pr...t.., group=row.names.x.)) +
    geom_histogram(bins = 20) +
    facet_wrap(~row.names.x.)

  fig$VolcanoPlot <- Model_Coeff %>%
    dplyr::filter(row.names.x. != "(Intercept)") %>%
    quantable::multigroupVolcano(
      effect = "Estimate",
      type = "Pr...t..",
      condition = "row.names.x.",
      label = hierarchyKey ,
      xintercept = c(-1, 1) ,
      colour = "isSingular" )

  forPairs <- Model_Coeff %>%
    dplyr::select(!!sym(hierarchyKey) , row.names.x. ,  Estimate ) %>%
    tidyr::spread(row.names.x.,Estimate )

  fig$Pairsplot_Coef <-  GGally::ggpairs(forPairs, columns=2:ncol(forPairs))

  fig$histogram_anova_p.values <- ggplot(data=modellingResult$Model_Anova, aes(x = Pr..F., group=rownames.x.)) +
    geom_histogram(bins = 20) +
    facet_wrap(~rownames.x.)


  return(fig)
}



#' Writes figures generated by `workflow_model_analyse_vis`
#'
#' used in p2901
#' @export
workflow_model_analyse_vis_write <- function(modelling_result,
                                             modelName,
                                             path,
                                             fig.width = 10 ,
                                             fig.height = 10){
  message("Writing figures into : ", path, "\n")

  pdf(file.path(path,paste0("Coef_Histogram_",modelName,".pdf")),
      width = fig.width, height = fig.height )
  print(modelling_result$histogram_coeff_p.values)
  dev.off()

  pdf(file.path(path,paste0("Coef_VolcanoPlot_",modelName,".pdf")),
      width = fig.width , height = fig.height)
  print(modelling_result$VolcanoPlot)
  dev.off()

  pdf(file.path(path, paste0("Coef_Pairsplot_",modelName,".pdf")),
      width = fig.width , height = fig.height)
  print(modelling_result$Pairsplot_Coef)
  dev.off()

  pdf(file.path(path,paste0("Anova_p.values_", modelName, ".pdf")),
      width = fig.width , height = fig.height)
  histogram_anova_p.values <-
    print(modelling_result$histogram_anova_p.values)
  dev.off()
}


#' p2621 workflow likelihood ratio test
#' @export
workflow_likelihood_ratio_test <- function(modelProteinF,
                                           modelName,
                                           modelProteinF_Int,
                                           modelName_Int,
                                           subject_Id = "protein_Id",
                                           path = NULL
){
  # Model Comparison
  reg <- inner_join(dplyr::select(modelProteinF, !!sym(subject_Id), starts_with("lmer_")),
                    dplyr::select(modelProteinF_Int, !!sym(subject_Id), starts_with("lmer_")) , by=subject_Id)

  reg <- reg %>% dplyr::mutate(modelComparisonLikelihoodRatioTest = map2(!!sym(paste0("lmer_", modelName)),
                                                                  !!sym(paste0("lmer_", modelName_Int)),
                                                                  .likelihood_ratio_test ))
  likelihood_ratio_test_result <- reg %>%
    dplyr::select(!!sym(subject_Id), modelComparisonLikelihoodRatioTest) %>% tidyr::unnest()
  likelihood_ratio_test_result <- likelihood_ratio_test_result %>%
    dplyr::rename(likelihood_ratio_test.pValue = modelComparisonLikelihoodRatioTest)


  if(!is.null(path)){
    fileName <- paste("hist_LRT_", modelName, "_", modelName_Int,".pdf", sep="")
    fileName <- file.path(path, fileName)
    message("writing figure : " , fileName , "\n")
    pdf(fileName)
    par(mfrow=c(2,1))
    hist(likelihood_ratio_test_result$likelihood_ratio_test.pValue, breaks=20)
    plot(ecdf(likelihood_ratio_test_result$likelihood_ratio_test.pValue))
    abline(v=c(0.01,0.05), col=c(3,2))
    dev.off()
  }

  return(likelihood_ratio_test_result)
}




