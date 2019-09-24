# mixed linear models ----

# Creating models from configuration ----


#' Create custom lmer model
#'
#' @export
#' @examples
#'
#' tmp <- make_custom_model_lmer("Intensity ~ condition + (1|peptide_Id)", model_name="random_example")
#' tmp$model_fun(get_formula=TRUE)
#' tmp$isSingular
#'
make_custom_model_lmer <- function( modelstr, model_name = "Model") {
  formula <- as.formula(modelstr)
  model_fun <- function(x, get_formula=FALSE){
    if(get_formula)
    {
      return(formula)
    }
    modelTest <- tryCatch(lmerTest::lmer( formula , data=x ),
                          error = function(e){print(e) ; return=NULL})
    return(modelTest)
  }
  res <- list(model_fun = model_fun,
              isSingular = lme4::isSingular,
              contrast_fun = my_contest,
              model_name = model_name,
              report_columns = c("p.value", "p.value.adjusted") )
  return(res)
}

#' Create custom ml model
#' @export
#' @examples
#' tmp <- make_custom_model_lm("Intensity ~ condition", model_name = "parallel design")
#' tmp$model_fun(get_formula=TRUE)
#' tmp$isSingular
make_custom_model_lm <- function( modelstr, model_name) {
  formula <- as.formula(modelstr)
  model_fun <- function(x, get_formula=FALSE){
    if(get_formula)
    {
      return(formula)
    }
    modelTest <- tryCatch(lm( formula , data=x ),
                          error = function(e){print(e) ; return=NULL})
    return(modelTest)
  }
  res <- list(model_fun = model_fun,
              isSingular = isSingular_lm,
              contrast_fun = my_contrast_V2,
              model_name = model_name,
              report_columns = c("moderated.p.value", "moderated.p.value.adjusted") )
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


# Fit the models to data ----



#' check if lm model is singular
#' @export
#'
isSingular_lm <- function(m){
  anyNA <- any(is.na(coefficients(m)))
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
#' library(tidyverse)
#' library(LFQService)
#' D <- LFQService::resultsV12954
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)",
#'   model_name = modelName)
#' pepIntensity <- D$pepIntensityNormalized
#' config <- D$config_pepIntensityNormalized
#' config$table$hkeysLevel()
#' mr <- model_analyse( pepIntensity,
#'  formula_randomPeptide,
#'  modelName,
#'  config$table$hkeysLevel())
#' get_complete_model_fit(mr$modelProtein)
model_analyse <- function(pepIntensity,
                          modelFunction,
                          modelName,
                          subject_Id = "protein_Id")
{
  pepIntensity %>%
    dplyr::group_by(!!!syms(subject_Id)) %>%
    tidyr::nest() -> nestProtein

  lmermodel <- "linear_model"

  nestProtein %>% dplyr::mutate(!!lmermodel := purrr::map(data, modelFunction$model_fun)) -> modelProtein

  modelProtein <- modelProtein %>% dplyr::mutate(!!"exists_lmer" := purrr::map_lgl(!!sym(lmermodel), function(x){!is.null(x)}))

  modelProteinF <- modelProtein %>% dplyr::filter( !!sym("exists_lmer") == TRUE)
  modelProteinF <- modelProteinF %>% dplyr::mutate(!!"isSingular" := purrr::map_lgl(!!sym(lmermodel), modelFunction$isSingular ))
  modelProteinF <- modelProteinF %>% dplyr::mutate(!!"df.residual" := purrr::map_dbl(!!sym(lmermodel), df.residual ))
  modelProteinF <- modelProteinF %>% dplyr::mutate(!!"sigma" := purrr::map_dbl( !!sym(lmermodel) , sigma))

  nrcoeff <- function(x){
    cc <- coefficients(x)
    if(class(cc) == "numeric"){
      return(length(cc))
    }else{
      return(ncol(cc[[1]]))
    }
  }

  nrcoeff_not_NA <- function(x){
    cc <- coefficients(x)
    if(class(cc) == "numeric"){
      return(sum(!is.na(cc)))
    }else{
      return(ncol(cc[[1]]))
    }
  }

  modelProteinF <- modelProteinF %>% dplyr::mutate(nrcoef = purrr::map_int(!!sym(lmermodel), nrcoeff))
  modelProteinF <- modelProteinF %>% dplyr::mutate(nrcoeff_not_NA = purrr::map_int(!!sym(lmermodel), nrcoeff_not_NA))

  #return(list(modelProtein = modelProtein, modelProteinF = modelProteinF))
  modelProteinF <- modelProteinF %>%
    dplyr::select_at(c(subject_Id,"isSingular", "df.residual","sigma" ,"nrcoef", "nrcoeff_not_NA") )
  modelProtein <- dplyr::left_join(modelProtein, modelProteinF)

  return(list(modelProtein = modelProtein,
              modelName = modelName
  ))
}



#' retrieve complete model.
#' @export
#'
get_complete_model_fit <- function(modelProteinF){
  modelProteinF <- modelProteinF %>% dplyr::filter(!!sym("exists_lmer") == TRUE)
  modelProteinF <- modelProteinF %>% dplyr::filter(nrcoeff_not_NA == max(nrcoeff_not_NA)) %>% dplyr::arrange(dplyr::desc(nrcoeff_not_NA))
  modelProteinF <- modelProteinF %>% dplyr::filter(df.residual > 0)
  return(modelProteinF)
}


#' summarize - compute anova and extract model coefficients from generated by `model_analyse`
#'
#' @export
#' @examples
#' rm(list=ls())
#' library(LFQService)
#' D <- LFQService::resultsV12954
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)",
#'    model_name=modelName)
#' pepIntensity <- D$pepIntensityNormalized
#' config <- D$config_pepIntensityNormalized
#' config$table$hkeysLevel()
#' modellingResult <- model_analyse( pepIntensity,
#'  formula_randomPeptide,
#'  modelName,
#'  config$table$hkeysLevel())
#' names(modellingResult)
#' tmp <- model_analyse_summarize(modellingResult$modelProtein, modelName)
#' names(tmp)
model_analyse_summarize <- function(modelProteinF,
                                    modelName = "Model",
                                    subject_Id = "protein_Id"){
  lmermodel <- "linear_model"

  modelProteinF <- get_complete_model_fit(modelProteinF)
  # modelProteinF <- modelProteinF %>% dplyr::filter(nrcoef == max(nrcoef))

  # Extract coefficients
  .coef_df <-  function(x){
    x <- coef(summary(x));
    x<- data.frame(row.names(x), x);
    return(x)
  }

  Model_Coeff <- modelProteinF %>%
    dplyr::mutate(!!"Coeffs_model" := purrr::map( !!sym(lmermodel),  .coef_df ))

  Model_Coeff <- Model_Coeff %>%
    dplyr::select(!!!syms(subject_Id), !!sym("Coeffs_model"), isSingular, nrcoef) %>%
    tidyr::unnest()

  # ANOVA
  .anova_df <- function(x){
    x <- anova(x)
    colnames(x) <- make.names(colnames(x))
    x <- data.frame(rownames(x), x)
    return(x)
  }

  Model_Anova <- modelProteinF %>% dplyr::mutate(!!"Anova_model" := purrr::map( !!sym(lmermodel),  .anova_df ))

  Model_Anova <- Model_Anova %>%
    dplyr::select(!!!syms(subject_Id), !!sym("Anova_model"), isSingular, nrcoef) %>%
    tidyr::unnest()


  return(list(
    modelName = modelName,
    Model_Coeff = Model_Coeff,
    fname_Model_Coeff =  paste0("Coef_",modelName, ".csv"),
    Model_Anova = Model_Anova,
    fname_Model_Anova =  paste0("ANOVA_",modelName,".csv" )
  ))
}


#' writes results of `model_analyse`, anova table and all the coefficients with parameters.
#' @export
model_analyse_summarize_write  <- function(modellingResult, path, all=FALSE){
  message("writing tables into :", path)
  if(all){
    lfq_write_table(modellingResult$Model_Coeff,
                    path = file.path( path,modellingResult$fname_Model_Coeff ))
  }
  lfq_write_table(modellingResult$Model_Anova,
                  path = file.path( path , modellingResult$fname_Model_Anova ))
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
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)",
#'    model_name = modelName)
#' modellingResult <-  model_analyse(
#'  D$pepIntensityNormalized,
#'  formula_randomPeptide,
#'  modelName,
#'  D$config_pepIntensityNormalized$table$hkeysLevel())
#' tmp <- model_analyse_summarize(modellingResult$modelProtein)
#' res <- model_analyse_summarize_vis(tmp,
#'  D$config_pepIntensityNormalized$table$hkeysLevel())
#'
model_analyse_summarize_vis <- function(modellingResult, subject_Id ="protein_Id") {
  Model_Coeff <- tidyr::unite(modellingResult$Model_Coeff, "subject_Id", subject_Id)
  Model_Anova <- tidyr::unite(modellingResult$Model_Anova, "subject_Id", subject_Id)
  modelName <- modellingResult$modelName
  fig <- list()

  ## Coef_Histogram
  fig$fname_histogram_coeff_p.values <- paste0("Coef_Histogram_",modelName,".pdf")
  fig$histogram_coeff_p.values <- ggplot(data = Model_Coeff, aes(x = Pr...t.., group=row.names.x.)) +
    geom_histogram(bins = 20) +
    facet_wrap(~row.names.x.)

  ## Coef_VolcanoPlot
  fig$fname_VolcanoPlot <- paste0("Coef_VolcanoPlot_",modelName,".pdf")
  fig$VolcanoPlot <- Model_Coeff %>%
    dplyr::filter(row.names.x. != "(Intercept)") %>%
    quantable::multigroupVolcano(
      effect = "Estimate",
      p.value = "Pr...t..",
      condition = "row.names.x.",
      label = "subject_Id" ,
      xintercept = c(-1, 1) ,
      colour = "isSingular" )

  ## Coef_Pairsplot
  forPairs <- Model_Coeff %>%
    dplyr::select(!!sym("subject_Id") , row.names.x. ,  Estimate ) %>%
    tidyr::spread(row.names.x.,Estimate )
  fig$fname_Pairsplot_Coef <- paste0("Coef_Pairsplot_",modelName,".pdf")
  fig$Pairsplot_Coef <-  GGally::ggpairs(forPairs, columns=2:ncol(forPairs))

  ## Anova_p.values
  fig$fname_histogram_anova_p.values <- paste0("Anova_p.values_", modelName, ".pdf")
  fig$histogram_anova_p.values <-  modellingResult$Model_Anova %>% dplyr::filter(rownames.x. != "Residuals") %>%
    ggplot( aes(x = Pr..F., group=rownames.x.)) +
    geom_histogram(bins = 20) +
    facet_wrap(~rownames.x.)

  return(fig)
}

#' Writes figures generated by `model_analyse_summarize_vis`
#'
#' used in p2901
#' @export
model_analyse_summarize_vis_write <- function(modelling_result,
                                              path,
                                              fig.width = 10 ,
                                              fig.height = 10,
                                              all = FALSE){
  if(all){
    fpath <- file.path(path, modelling_result$fname_histogram_coeff_p.values)
    message("Writing figure into : ", fpath, "\n")
    pdf(fpath, width = fig.width, height = fig.height )
    print(modelling_result$histogram_coeff_p.values)
    dev.off()

    fpath <- file.path(path, modelling_result$fname_VolcanoPlot)
    message("Writing figure into : ", fpath, "\n")
    pdf(fpath,
        width = fig.width , height = fig.height)
    print(modelling_result$VolcanoPlot)
    dev.off()

    fpath <- file.path(path, modelling_result$fname_Pairsplot_Coef)
    message("Writing figure into : ", fpath, "\n")
    pdf(fpath, width = fig.width , height = fig.height)
    print(modelling_result$Pairsplot_Coef)
    dev.off()
  }
  fpath <- file.path(path, modelling_result$fname_histogram_anova_p.values)
  message("Writing figure into : ", fpath, "\n")
  pdf(fpath, width = fig.width , height = fig.height)
  print(modelling_result$histogram_anova_p.values)
  dev.off()
}

#' workflow_model_analyse
#' @export
#' @references function with paramter path
#' @examples
#' D <- LFQService::resultsV12954
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ Condition + (1 | peptide_Id)",
#'   model_name = modelName)
#' modellingResult <-  workflow_model_analyse(D$pepIntensityNormalized,
#'  formula_randomPeptide,
#'   modelName,
#'  subject_Id = D$config_pepIntensityNormalized$table$hkeysLevel())
#' reslist <- modellingResult()
workflow_model_analyse <- function(data,
                                   modelFunction,
                                   modelName = "Model",
                                   subject_Id = "protein_Id"){

  modellingResult <- model_analyse(data,
                                   modelFunction,
                                   modelName = modelName,
                                   subject_Id)

  # delay write
  res_fun <- function(path = NULL, all=FALSE){
    summaryResult <- model_analyse_summarize(modellingResult$modelProtein,
                                             modelName = modelName,
                                             subject_Id = subject_Id)
    visualization <- model_analyse_summarize_vis(summaryResult, subject_Id)

    if(!is.null(path)){
      model_analyse_summarize_write(summaryResult, path, all=all)
      model_analyse_summarize_vis_write(visualization, path,all=all)
    }
    return(list(modellingResult = modellingResult,
                summaryResult = summaryResult,
                visualization = visualization))
  }
  return(res_fun)
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
  reg <- dplyr::inner_join(dplyr::select(modelProteinF, !!sym(subject_Id), "linear_model"),
                           dplyr::select(modelProteinF_Int, !!sym(subject_Id), "linear_model") , by=subject_Id)

  reg <- reg %>% dplyr::mutate(modelComparisonLikelihoodRatioTest = map2(!!sym("linear_model.x"),
                                                                         !!sym("linear_model.y"),
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
#' library(LFQService)
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
  head(data)
  head(ran)
  ran <- dplyr::inner_join(data, ran, by=randeffect)

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
    ran <- dplyr::inner_join(data, ran, by=rand_i)
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
  #ran_res <- dplyr::inner_join(ran_res, ran, by=rand_i)
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
#' m <- LFQService::interactionModel_p1807
#' plot_lmer_predicted_interactions(plot_lmer_model_and_data(m,"dumm"),m)
plot_lmer_predicted_interactions <- function(gg, m){
  cm <- .lmer4_coeff_matrix(m)
  xstart_end <- data.frame(xstart = rownames(cm$mm), xend = rownames(cm$mm))
  ystart_end <- data.frame(xend = rownames(cm$mm), ystart =rep(0, nrow(cm$mm)),
                           yend = cm$mm %*% cm$coeffs)
  segments <- dplyr::inner_join(xstart_end, ystart_end, by="xend")
  gg <- gg + geom_segment(aes(x = xstart, y = ystart , xend = xend, yend =yend), data=segments, color = "blue", arrow=arrow())
  return(gg)
}

#' Make model plot with title - protein Name.
#' @export
#' @examples
#' m <- LFQService::interactionModel_p1807
#' plot_lmer_model_and_data(m,"dumm")
#'
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
#' all.equal(linfct$linfct_factors["CelltypeCMP/MEP",] ,
#'  apply(linfct$linfct_interactions[grep("CelltypeCMP/MEP", rownames(linfct$linfct_interactions)),],2, mean))
#' linfct$linfct_interactions
#' #}
#'
#' m <- lm(Petal.Width ~ Species, data=iris)
#' linfct_from_model(m)
linfct_from_model <- function(m, as_list = TRUE){

  cm <- .lmer4_coeff_matrix(m)
  cm_mm <- cm$mm[order(rownames(cm$mm)),]

  dd <- .coeff_weights_factor_levels(cm_mm)
  dd_m <- dd %>% dplyr::select(-factor_level) %>% data.matrix()
  rownames(dd_m) <- dd$factor_level
  dd_m <- dd_m[order(rownames(dd_m)),]
  res <- list(linfct_factors = dd_m , linfct_interactions = cm_mm)

  if(as_list){
    return(res)
  }else{
    do.call( rbind, res)
  }
}

#' linfct_matrix_contrasts
#' @export
#' @examples
#' m <- LFQService::basicModel_p1807
#' linfct <- linfct_from_model(m,as_list=FALSE)
#' linfct
#'
#' Contrasts <- c("CMP/MEP - HSC" = "`CelltypeCMP/MEP` - `CelltypeHSC`",
#' "someWeird" = "`class_therapyc.NO:CelltypeCMP/MEP` - `class_therapyp.HU:CelltypeCMP/MEP`")
#' linfct_matrix_contrasts(linfct, Contrasts )
linfct_matrix_contrasts<- function(linfct , contrasts){
  linfct <- t(linfct)
  df <- as_tibble(linfct, rownames = "interaction")
  make_contrasts <- function(data,
                        contrasts)
  {
    cnams <- setdiff(colnames(data),"interaction")
    for(i in 1:length(contrasts)){
      message(names(contrasts)[i], "=", contrasts[i],"\n")
      data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
    }

    res <- data %>% dplyr::select(-one_of(cnams))
    return(res)
  }

  res <- make_contrasts(df, contrasts )
  res <- column_to_rownames(res,"interaction")
  res <- t(res)
  return(res)
}


#' create all possible contrasts
#' @examples
#' m <- LFQService::basicModel_p1807
#' m
#' linfct <- linfct_from_model(m)
#'
#' xl <- LFQService:::.linfct_all_possible_contrasts(linfct$linfct_factors)
#' xx <- LFQService:::.linfct_all_possible_contrasts(linfct$linfct_interactions)
#'
.linfct_all_possible_contrasts <- function( lin_int ){
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
#'
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
  ffac <- ffac[!grepl(":",ffac)] # remove interactions
  linfct_factors <- linfct_from_model(m)$linfct_factors

  factorLevels <- rownames(linfct_factors)
  res <- vector(length(ffac), mode = "list")
  for(i in 1:length(ffac)){
    fac <- ffac[i]
    idx <- grep(fac, factorLevels)
    linfct_m <- linfct_factors[idx,]
    res[[i]] <- .linfct_all_possible_contrasts(linfct_m)
  }
  res <- do.call(rbind, res)
  return(res)
}

# Computing contrasts helpers -----

#' apply multcomp::glht method to linfct
#'
#' @export
#' @examples
#'
#' mb <- LFQService::basicModel_p1807
#' linfct <- linfct_from_model(mb)
#' names(linfct)
#' my_glht(mb, linfct$linfct_factors)
#'
#' m <- LFQService::modellingResult_A$modelProtein$linear_model[[1]]
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_glht(m, linfct)
#'
my_glht <- function(model, linfct , sep=TRUE ) {
  if(!class(model) == "lm") # fixes issue of mutlcomp not working on factors of class character
  {
    warning("USE ONLY WITH LM models ", class(model))
    if(length(lme4::fixef(model)) != ncol(linfct) ){
      return(NA) # catch rank defficient
    }
  }else{
    if(isSingular_lm(model)){
      return(NA)
    }
    model$model <- as.data.frame(unclass(model$model))
  }
  if(sep){
    res <- list()
    for(i in 1:nrow(linfct)){
      x <- multcomp::glht(model, linfct=linfct[i,,drop=FALSE])
      RHS <- broom::tidy(confint(x)) %>% dplyr::select(-estimate)

      RHS$df <- x$df
      RHS$sigma <- sigma(model)

      x <- dplyr::inner_join(broom::tidy(summary(x)),RHS,by = c("lhs", "rhs")) %>% dplyr::select(-rhs)
      res[[i]] <- x
    }
    res <- bind_rows(res)
    return(res)
  }else{
    x <- multcomp::glht(model, linfct = linfct)
    RHS <- broom::tidy(confint(x)) %>% dplyr::select(-estimate)
    RHS$df <- x$df
    RHS$sigma <- sigma(model)
    res <- dplyr::inner_join(broom::tidy(summary(x)), RHS, by = c("lhs", "rhs")) %>%
      dplyr::select(-rhs)
    return(res)
  }
}

#' compute contrasts for full models
#'
#' @export
#' @examples
#'
#' m <- LFQService::modellingResult_A$modelProtein$linear_model[[1]]
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_glht(m, linfct)
#' my_contrast(m, linfct, confint = 0.95)
#' my_contrast(m, linfct, confint = 0.99)
#'
my_contrast <- function(m,
                        linfct,
                        coef = coefficients(m),
                        Sigma.hat = vcov(m), confint = 0.95){

  df <- df.residual(m)
  sigma <- sigma(m)

  estimate <- linfct %*% t(t(coef))

  if(df > 0){
    std.error <- sqrt(diag(linfct %*% Sigma.hat %*% t(linfct)))
    statistic <- estimate / std.error

    #p.value <- pt(-abs(statistic), df = df) * 2

    p.value <- pt(abs(statistic), df = df, lower.tail = FALSE) * 2
    prqt <- -qt((1-confint)/2, df=df)
    conf.low <- estimate  - prqt * std.error
    conf.high <- estimate + prqt * std.error

  }else{
    std.error <- NA
    statistic <- NA
    p.value <- NA
    conf.low <- NA
    conf.high <- NA
  }

  res <- data.frame(lhs = rownames(linfct),
                    sigma = sigma,
                    df = df,
                    estimate = estimate,
                    std.error = std.error,
                    statistic = statistic ,
                    p.value = p.value,
                    conf.low= conf.low,
                    conf.high =conf.high, stringsAsFactors = FALSE)
  return(res)
}

#' handles incomplete models by setting coefficients to 0
#' @export
#' @examples
#' m <- LFQService::modellingResult_A$modelProtein$linear_model[[1]]
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_glht(m, linfct)
#' my_contrast_V1(m, linfct, confint = 0.95)
#' my_contrast_V1(m, linfct, confint = 0.99)
my_contrast_V1 <- function(incomplete, linfct,confint = 0.95){
  Sigma.hat <- vcov(incomplete)
  Sigma.hat[is.na(Sigma.hat)] <- 0
  coef <- coefficients(incomplete)
  coef[is.na(coef)] <- 0
  res <- my_contrast(incomplete, linfct,
                     coef = coef,
                     Sigma.hat = Sigma.hat,
                     confint = confint)
  return(res)
}

#' handles incomplete models
#'
#' only keeps non NA coefficients.
#'
#' @export
#' @examples
#' m <- LFQService::modellingResult_A$modelProtein$linear_model[[1]]
#' linfct <- linfct_from_model(m)$linfct_factors
#' my_contrast_V2(m, linfct, confint = 0.95)
#' my_contrast_V2(m, linfct, confint = 0.99)
#'
my_contrast_V2 <- function(m, linfct,confint = 0.95){
  Sigma.hat <- vcov(m)
  coef <- coefficients(m)
  res <- vector(nrow(linfct), mode="list")
  for(i in 1:nrow(linfct)){
    linfct_v <- linfct[i,,drop=FALSE]
    idx <- which(linfct_v != 0)
    nam <- colnames(linfct_v)[idx]

    if(all(nam %in% names(coef))){
      linfct_v_red <- linfct_v[, nam, drop=FALSE]
      Sigma.hat_red <- Sigma.hat[nam,nam,drop=FALSE]
      coef_red <- coef[nam]
      stopifnot(all.equal(colnames(linfct_v_red),colnames(Sigma.hat_red)))
      stopifnot(all.equal(colnames(linfct_v_red),names(coef_red)))
      res[[i]] <- my_contrast(m,linfct_v_red,
                              coef=coef_red,
                              Sigma.hat = Sigma.hat_red,confint = confint)
    }else{
      res[[i]] <-  data.frame(lhs = rownames(linfct_v),
                              sigma = sigma(m),
                              df = df.residual(m),
                              estimate = NA,
                              std.error = NA,
                              statistic = NA ,
                              p.value = NA,
                              conf.low= NA,
                              conf.high =NA,
                              stringsAsFactors = FALSE)
    }
  }
  res <- bind_rows(res)
  return(res)
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
#' #my_glht(mb, linfct$linfct_factors)
#' #my_glht(mb, linfct$linfct_interactions)
my_contest <- function(model, linfct){
  if(length(lme4::fixef(model)) != ncol(linfct) ){
    warning("Model is rank deficient!")
    return(NA) # catch rank defficient
  }
  res <- lmerTest::contest(model, linfct, joint = FALSE, confint = TRUE)
  res <- as_tibble(res, rownames="lhs")
  res$sigma <- sigma(model)
  res <- res %>% rename(estimate = Estimate,
                        std.error = "Std. Error",
                        statistic="t value",
                        p.value = "Pr(>|t|)",
                        conf.low = "lower",
                        conf.high = "upper")
  return(res)
}


# computing contrast ----

#' pivot model contrasts matrix to wide format produced by `contrasts_linfct` and ...
#' @export
#' @examples
#' dd <- LFQService::factor_levelContrasts
#' head(dd)
#' tmp <- pivot_model_contrasts_2_Wide(dd, subject_Id = "Compound")
#' colnames(tmp)
pivot_model_contrasts_2_Wide <- function(modelWithInteractionsContrasts,
                                         subject_Id = "protein_Id",
                                         columns = c("estimate", "p.value","p.value.adjusted")){

  m_spread <- function(longContrasts, subject_Id, column ){
    longContrasts %>%
      dplyr::select_at(c(subject_Id, "isSingular", "lhs",  column)) %>%
      dplyr::mutate(lhs = glue::glue('{column}.{lhs}')) %>%
      tidyr::spread(lhs, !!sym(column) ) -> res
    return(res)
  }
  res <- list()
  for(column in columns){
    res[[column]] <- m_spread(modelWithInteractionsContrasts,subject_Id,column)
  }
  res <- res %>% reduce(left_join, by = c(subject_Id,"isSingular"))
  return(res)
}
#' compute group averages
#'
#' used in p2621, p2109
#'
#' @export
#' @examples
#'
#' modelSummary_A <- LFQService::modellingResult_A
#' m <- get_complete_model_fit(modelSummary_A$modelProtein)
#'
#'
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#' factor_contrasts
#'
#' factor_levelContrasts <- contrasts_linfct( m,
#'         factor_contrasts,
#'         subject_Id = "Compound",
#'         contrastfun = LFQService::my_contrast_V2)
#'
#' #usethis::use_data(factor_levelContrasts, overwrite=TRUE)
#'
#' models_interaction <- LFQService::models_interaction
#' #models_interaction$modelProtein <- dplyr::rename(models_interaction$modelProtein, linear_model = lmer_f_Condition_r_peptid_r_patient)
#' #usethis::use_data(models_interaction,overwrite = TRUE)
#' m <- get_complete_model_fit(models_interaction$modelProtein)
#' factor_contrasts <- linfct_factors_contrasts( m$linear_model[[1]])
#' m
#' factor_levelContrasts <- contrasts_linfct( m,
#'                            factor_contrasts,
#'                        subject_Id = "protein_Id")
#' head(factor_levelContrasts)
#' plot(factor_levelContrasts$df, factor_levelContrasts$df.residual.model )
#' plot(factor_levelContrasts$df.residual.model , factor_levelContrasts$df - factor_levelContrasts$df.residual.model )
#'
contrasts_linfct <- function(models,
                             linfct,
                             subject_Id = "protein_Id" ,
                             contrastfun = LFQService::my_contest){
  #computeGroupAverages
  modelcol <- "linear_model"
  models <- models %>% dplyr::filter(exists_lmer == TRUE)


  interaction_model_matrix <- models %>%
    dplyr::mutate(contrast = map(!!sym(modelcol) , contrastfun , linfct = linfct ))


  mclass <- function(x){
    class(x)[1]
  }

  interaction_model_matrix %>%
    dplyr::mutate(classC = map_chr(contrast,mclass)) %>%
    dplyr::filter(classC != "logical") -> interaction_model_matrix

  contrasts <- interaction_model_matrix %>%
    dplyr::select_at( c(subject_Id, "contrast") ) %>% tidyr::unnest()

  modelInfos <- models %>%
    dplyr::select_at(c(subject_Id, "isSingular",
                       "sigma.model" = "sigma",
                       "df.residual.model" = "df.residual" )) %>% distinct()
  contrasts <- dplyr::inner_join(contrasts, modelInfos, by=subject_Id)

  # adjust
  contrasts <- contrasts %>% group_by_at("lhs") %>%
    dplyr::mutate(p.value.adjusted = p.adjust(p.value, method="BH")) %>% ungroup()

  return(contrasts)
}



.multigroupVolcano <- function (data,
                                effect = "fc",
                                p.value = "p.adjust",
                                condition = "condition",
                                colour = "colour",
                                xintercept = c(-2,2),
                                pvalue = 0.05,
                                text = NULL,
                                ablines = data.frame(fc = c(0, 0), p = c(0.01, 0.05), Area = c("p=0.01", "p=0.05")),
                                scales = "fixed",
                                maxNrOfSignificantText = 20)
{
  colname = paste("-log10(", p.value, ")", sep = "")
  p <- ggplot(data, aes_string(x = effect, y = colname, color = colour, text = "label")) +
    geom_point(alpha = 0.5)
  p <- p + scale_colour_manual(values = c("black", "green",
                                          "blue", "red"))
  p <- p + facet_wrap(as.formula(paste("~", condition)),
                      scales = scales) + labs(y = colname)
  ablines$neg_log10p <- -log10(ablines$p)
  p <- p + geom_abline(data = ablines, aes_string(slope = "fc",
                                                  intercept = "neg_log10p", colour = "Area")) +
    geom_vline(xintercept = xintercept, linetype = "dashed",
               colour = "red") +
    theme_light()

  return(p)
}



#' visualize output of `contrasts_linfct``
#' @export
#'
contrasts_linfct_vis <- function(contrasts,
                                 modelName = "Model",
                                 prefix = "Contrasts",
                                 subject_Id = "protein_Id",
                                 columns = c("p.value","p.value.adjusted"),
                                 fc = 1){
  res <- list()
  contrasts %>% tidyr::unite("label", subject_Id, sep="~", remove=FALSE) -> contrasts
  # add histogram of p-values
  for(column in columns){
    fig <- list()
    name <- paste0(prefix,"_Histogram_",column)
    fig$fname <- paste0(name, "_", modelName )
    fig$fig <- ggplot(data=contrasts, aes(x = !!sym(column))) +
      geom_histogram(bins = 20) +
      facet_wrap(~lhs)
    res[[name]] <- fig
  }
  # add volcano plots
  for(column in columns){
    fig <- list()
    name <- paste0(prefix,"_Volcano_",column)
    fig$fname <- paste0(name, "_", modelName )


    fig$fig <- LFQService:::.multigroupVolcano(contrasts,
                                  effect = "estimate",
                                  p.value = column,
                                  condition = "lhs",
                                  text = "label",
                                  xintercept = c(-fc, fc),
                                  colour = "isSingular",
                                  scales="free_y")

    fig$plotly <- contrasts %>% plotly::highlight_key(~label) %>%
      LFQService:::.multigroupVolcano(.,
                         effect = "estimate",
                         p.value = "p.value",
                         condition = "lhs",
                         text = "label",
                         xintercept = c(-fc, fc),
                         colour = "isSingular",
                         scales="free_y") %>%
      plotly::ggplotly(tooltip = "label")

    res[[name]] <- fig
  }

  # add histogram of fold changes
  {
    fig <- list()
    name <- paste0(prefix,"_Histogram_FC_esimate")
    fig$fname <- paste0(name, "_", modelName )
    fig$fig <- ggplot(data=contrasts, aes(x = !!sym("estimate"))) +
      geom_histogram(breaks = seq(floor(min(contrasts$estimate)),ceiling(max(contrasts$estimate)), by=0.5)) +
      facet_wrap(~lhs)
    res[[name]] <- fig
  }
  # MA plot
  {
    ma_plot <- function(x, fc = 1){
      x <- ggplot(x , aes(x = (c1+c2)/2, y = estimate, text = !!sym("label"), colour = !!sym("isSingular"))) +
        geom_point(alpha = 0.5) + scale_colour_manual(values = c("black", "red")) +
        facet_wrap(~lhs) + theme_light() +
        geom_hline(yintercept = c(-fc, fc), linetype = "dashed",colour = "red")
      return(x)
    }

    if(!is.null(contrasts$c1) && !is.null(contrasts$c2)){
      fig <- list()
      name <- paste0(prefix,"_MA_FC_estimate")
      fig$fname <- paste0(name, "_", modelName )
      # pdf version
      fig$fig <- contrasts %>% ma_plot(fc = fc)

      # html version
      fig$plotly  <- contrasts %>%
        plotly::highlight_key(~label) %>%
        ma_plot() %>%
        plotly::ggplotly(tooltip = "label")

      res[[name]] <- fig
    }
  }
  return(res)
}

#' helper function to write the result of `contrasts_linfct_vis`
#'
#' used in p2901
#'
#' @export
contrasts_linfct_vis_write <- function(fig_list,
                                       path,
                                       fig.width = 10,
                                       fig.height = 10,
                                       format = c("pdf","html")){
  format <- match.arg(format)
  if(!is.null(path)){
    for(fig in fig_list){

      fpath <- file.path(path,paste0(fig$fname,".", format))


      if(format == "pdf"){
        message("Writing: ",fpath,"\n")
        pdf(fpath, width = fig.width, height = fig.height)
        print(fig$fig)
        dev.off()
      }else if(format == "html"){
        if(!is.null(fig$plotly)){
          message("Writing: ",fpath,"\n")
          htmlwidgets::saveWidget(widget=fig$plotly, fig$fname, selfcontained = TRUE)
          file.rename(fig$fname, fpath)
        }
      }
    }
  }
}

#' Do contrast
#' @export
#' @examples
#'
workflow_contrasts_linfct <- function(models,
                                      contrasts,
                                      config,
                                      modelName = "Model",
                                      prefix = "Contrasts",
                                      contrastfun = LFQService::my_contest )
{
  if(class(contrasts) == "matrix"){
    linfct_A <- contrasts
  }else{
    models <- models %>% dplyr::filter(exists_lmer == TRUE)
    m <- get_complete_model_fit(models)
    linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
    linfct_A <- linfct_matrix_contrasts(linfct, contrasts)
  }

  subject_Id <- config$table$hkeysLevel()
  contrast_result <- contrasts_linfct(models,
                                      linfct_A,
                                      subject_Id = subject_Id,
                                      contrastfun = contrastfun )

  contrast_result <- moderated_p_limma_long(contrast_result)
  subject_Id <- subject_Id
  prefix <- prefix
  modelName <- modelName

  res_fun <- function(path = NULL, columns = c("p.value",
                                               "p.value.adjusted",
                                               "moderated.p.value",
                                               "moderated.p.value.adjusted")){

    visualization <- contrasts_linfct_vis(contrast_result,
                                          modelName ,
                                          prefix = prefix,
                                          subject_Id = subject_Id,
                                          columns = columns
    )

    relevant_columns <- c("lhs", "sigma", "df", "isSingular", "estimate", "conf.low", "conf.high") # other relevant columns.
    contrast_minimal <- contrast_result %>% dplyr::select(subject_Id, relevant_columns, columns )

    contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                   subject_Id = subject_Id,
                                                   columns=c("estimate", columns))

    if(!is.null(path)){
      if(FALSE){
        contrasts_linfct_write(contrast_minimal,
                               config,
                               path=path,
                               modelName = modelName,
                               prefix = prefix,
                               columns = c("estimate", columns))
      }

      contrasts_linfct_vis_write(visualization, path=path)
    }

    res <- list(contrast_result = contrast_result,
                contrast_minimal = contrast_minimal,
                contrasts_wide = contrasts_wide,
                visualization = visualization,
                modelName = modelName,
                prefix = prefix)

    invisible(res)
  }
  return( res_fun )
}


#' Do contrast
#' @export
#' @examples
#'
workflow_contrasts_linfct_V2 <- function(models,
                                         contrasts,
                                         config,
                                         modelName = "Model",
                                         prefix = "Contrasts",
                                         contrastfun = LFQService::my_contest )
{

  # extract contrast sides
  tt <- contrasts[grep("-",contrasts)]
  tt <- tibble(lhs = names(tt) , contrast= tt)
  tt <- tt %>% mutate(contrast = gsub("[` ]","",contrast)) %>%
    tidyr::separate(contrast, c("c1", "c2"), sep="-")



  models <- models %>% dplyr::filter(exists_lmer == TRUE)
  m <- get_complete_model_fit(models)
  linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
  linfct_A <- linfct_matrix_contrasts(linfct, contrasts)


  subject_Id <- config$table$hkeysLevel()
  contrast_result <- contrasts_linfct(models,
                                      rbind(linfct, linfct_A),
                                      subject_Id = subject_Id,
                                      contrastfun = contrastfun )

  xx <- contrast_result %>% dplyr::select(subject_Id, "lhs", "estimate")
  xx <- xx %>% pivot_wider(names_from = "lhs", values_from = "estimate")

  contrast_result <- contrast_result %>% dplyr::filter(lhs %in% names(contrasts))


  get_contrast_cols <- function(i, contrast_results , contrast_table , subject_ID ){
    data.frame(lhs = contrast_table[i, "lhs"],
               dplyr::select_at(contrast_results, c( subject_ID ,unlist(contrast_table[i,c("c1","c2")]))),
               c1_name = contrast_table[i,"c1", drop=T],
               c2_name = contrast_table[i,"c2", drop=T], stringsAsFactors = FALSE)
  }

  contrast_sides <- purrr::map_df(1:nrow(tt), get_contrast_cols, xx, tt, subject_Id)

  contrast_result <- inner_join(contrast_sides,contrast_result)


  contrast_result <- moderated_p_limma_long(contrast_result)
  subject_Id <- subject_Id
  prefix <- prefix
  modelName <- modelName

  res_fun <- function(path = NULL, columns = c("p.value",
                                               "p.value.adjusted",
                                               "moderated.p.value",
                                               "moderated.p.value.adjusted")){

    visualization <- contrasts_linfct_vis(contrast_result,
                                          modelName ,
                                          prefix = prefix,
                                          subject_Id = subject_Id,
                                          columns = columns
    )

    relevant_columns <- c("lhs", "c1_name", "c1", "c2_name", "c2", "sigma", "df", "isSingular", "estimate", "conf.low", "conf.high") # other relevant columns.
    contrast_minimal <- contrast_result %>% dplyr::select(subject_Id, relevant_columns, columns )

    contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                   subject_Id = subject_Id,
                                                   columns=c("estimate", columns))

    if(!is.null(path)){
      if(FALSE){
        contrasts_linfct_write(contrast_minimal,
                               config,
                               path=path,
                               modelName = modelName,
                               prefix = prefix,
                               columns = c("estimate", columns))
      }

      contrasts_linfct_vis_write(visualization, path=path, format = "pdf")
      contrasts_linfct_vis_write(visualization, path=path, format = "html")
    }

    res <- list(contrast_result = contrast_result,
                contrast_minimal = contrast_minimal,
                contrasts_wide = contrasts_wide,
                visualization = visualization,
                modelName = modelName,
                prefix = prefix)

    invisible(res)
  }
  return( res_fun )
}

# LIMMA ----

#' Moderate p-values - limma approach
#' @export
moderated_p_limma <- function(mm, df = "df"){
  sv <- limma::squeezeVar(mm$sigma^2, df=mm[[df]])
  sv <- as_tibble(sv)
  sv <- sv %>% setNames(paste0('moderated.', names(.)))
  mm <- bind_cols(mm, sv)
  mm <- mm %>% dplyr::mutate(moderated.statistic  =  statistic * sigma /  sqrt(moderated.var.post))
  mm <- mm %>% dplyr::mutate(moderated.df.total = !!sym(df) + moderated.df.prior)
  mm <- mm %>% dplyr::mutate(moderated.p.value = 2*pt( abs(moderated.statistic), df=moderated.df.total, lower.tail=FALSE) )
  mm <- mm %>% dplyr::mutate(moderated.p.value.adjusted = p.adjust(moderated.p.value, method="BH")) %>% ungroup()
  return(mm)
}

#' Moderate p-value for long table
#' @param mm result of `contrasts_linfct``
#' @param group_by_col colnames with contrast description - default 'lhs'
#' @export
#' @examples
#' library(LFQService)
#' modelSummary_A <- LFQService::modellingResult_A
#' m <- get_complete_model_fit(modelSummary_A$modelProtein)
#' factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
#' factor_levelContrasts <- contrasts_linfct( modelSummary_A$modelProtein,
#'                                                    factor_contrasts,
#'                                                    subject_Id = "Compound",
#'                                                    contrastfun = my_contrast_V2)
#'
#' mmm <- moderated_p_limma_long(factor_levelContrasts, group_by_col = "lhs")
#' plot(mmm$p.value, mmm$moderated.p.value, log="xy")
#' abline(0,1, col=2)
#'
#' # updating lmer model
#' models_interaction <- LFQService::models_interaction
#'
#' m <- get_complete_model_fit(models_interaction$modelProtein)
#' factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
#'
#' factor_levelContrasts <- contrasts_linfct(m,
#'                                          factor_contrasts,
#'                                          subject_Id = "protein_Id")
#'
#' mmm <- moderated_p_limma_long(factor_levelContrasts, group_by_col = "lhs")
#' head(mmm)
#' plot(mmm$p.value, mmm$moderated.p.value, log="xy")
#' abline(0,1, col=2)
#'
moderated_p_limma_long <- function( mm , group_by_col = "lhs"){
  dfg <- mm %>% group_by_at(group_by_col) %>% group_split()
  xx <- purrr::map_df(dfg, moderated_p_limma)
  return(xx)
}

#' write results of `contrasts_linfct`
#' @export
#'
contrasts_linfct_write <- function(results,
                                   config,
                                   path,
                                   modelName = "Model",
                                   prefix = "Contrasts",
                                   columns = c("estimate", "p.value", "p.value.adjusted")){

  subject_Id <- config$table$hkeysLevel()

  if(!is.null(path)){
    fileLong <- file.path(path,paste0(prefix, "_", modelName, ".csv"))
    message("Writing: ", fileLong, "\n")
    lfq_write_table(separate_hierarchy(results, config) , path = fileLong)
    fileWide <- file.path(path,paste0(prefix, "_", modelName, "_PIVOT.csv"))
    message("Writing: ", fileWide, "\n")
    resultswide <- pivot_model_contrasts_2_Wide(results,
                                                subject_Id = subject_Id,
                                                columns=columns)
    lfq_write_table(separate_hierarchy(resultswide, config), path = fileWide)
  }
}

# HELPER ----

#' get coefficients from all models
#' @export
#'
get_model_coefficients <- function(modeldata, config){
  l_coeff <- function(m){
    if(!is.null(m)){
      res <- as.numeric(coefficients(m))
      return(res)
    }
    return(numeric())
  }

  n_coeff <- function(m){
    if(!is.null(m)){
      res <- names(coefficients(m))
      return(res)
    }
    return(character())
  }

  xx <- modeldata %>%
    dplyr::mutate(coefficients_values = purrr::map(linear_model, l_coeff)) %>%
    dplyr::mutate(coefficients_names = purrr::map(linear_model, n_coeff))

  xxs <- xx %>% dplyr::select( config$table$hkeysLevel(),
                               coefficients_values,
                               coefficients_names)
  xxxn<-xxs %>% unnest()
  xxcoef <- xxxn %>% spread(coefficients_names,coefficients_values)
  return(xxcoef)
}


# ROPECA ----

#' p-value of protein from p.value of the median fold change peptide.
#' @param max.n limit number of peptides per protein.
#' @export
#' @examples
#' plot(get_p_values_pbeta(0.1,1:10))
#' abline(h=.05,col=2)
#' plot(get_p_values_pbeta(0.3,1:30))
#' abline(h=.05,col=2)
#' plot(get_p_values_pbeta(rep(0.1,30),rep(3,30)))
#'
get_p_values_pbeta <- function(median.p.value, n.obs , max.n = 10){
  n.obs <- pmin(n.obs, max.n)

  shape1 <- (n.obs/2 + 0.5)
  shape2 <- (n.obs - (n.obs/2 + 0.5) + 1)

  stopifnot(shape1 == shape2)
  res.p.value <- pbeta(median.p.value,
                       shape1 = shape1,
                       shape2 = shape2)
  return(res.p.value)
}

# TODO delete
.getMedianIDX <- function(nrows){
  if(nrows%%2 == 0){
    idx <- c(nrows/2 , nrows/2 +1)
  }else{
    idx <- ceiling(nrows/2)
  }
  return(idx)
}
# TODO delete
.summarize_y_by_x_median <- function(x,y){
  not.na.idx <- which(!is.na(x))
  x <- x[not.na.idx]
  y <- y[not.na.idx]
  ord <- order(x)
  y <- y[ord]
  idx <- .getMedianIDX(length(y))
  res <- mean(y[idx])
  return(res)
}


#' compute protein level fold changes and p.values (using beta distribution)
#' takes p-value of the scaled p-value
#'
#' @export
#'
#' @examples
#'
#' library(LFQService)
#' nrPep <- 10000
#' nrProtein <- 800
#' p.value <- runif(nrPep)
#' estimate <- sample(c(-1,1),nrPep, replace = TRUE)
#' protein_Id <- sample(1:800, size = nrPep,
#'   replace=TRUE, prob = dexp(seq(0,5,length=800)))
#'
#' plot(table(table(protein_Id)))
#'
#' testdata <- data.frame(lhs = "contrast1", protein_Id = protein_Id,
#'   estimate = estimate, p.value = p.value )
#' xx <- summary_ROPECA_median_p.scaled(testdata,
#'
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value="p.value",
#'                                     max.n = 10)
#' colnames(xx)
#' hist(testdata$p.value)
#' hist(xx$median.p.scaled, breaks=20)
#' hist(xx$median.p, breaks=20)
#' hist(xx$beta.based.significance, breaks = 20)
#'
summary_ROPECA_median_p.scaled <- function(contrasts_data,
                                           contrast = "lhs",
                                           subject_Id = "protein_Id",
                                           estimate = "estimate",
                                           p.value="moderated.p.value",
                                           max.n = 10){
  addscaled.p <- contrasts_data %>%
    dplyr::mutate(scaled.p = ifelse(!!sym(estimate) > 0, 1-!!sym(p.value) , !!sym(p.value)-1))

  summarized.protein <- addscaled.p %>% group_by_at(c(subject_Id, contrast))

  summarized.protein <- summarized.protein %>%
    summarize(
      n=n(),
      n_not_na = sum(!is.na(!!sym(estimate))),
      median.estimate = median(!!sym(estimate), na.rm=TRUE),
      sd.estimate = sd(!!sym(estimate), na.rm=TRUE),
      median.p.scaled = median(scaled.p, na.rm=TRUE))

  # scale it back here.
  summarized.protein <- summarized.protein %>%
    dplyr::mutate(median.p = 1- abs(median.p.scaled))
  summarized.protein <- summarized.protein %>%
    dplyr::mutate(beta.based.significance = get_p_values_pbeta(median.p  , n_not_na, max.n = max.n))
  summarized.protein <- summarized.protein %>%
    dplyr::mutate(n.beta = pmin(n_not_na, max.n))
  return(summarized.protein)
}

