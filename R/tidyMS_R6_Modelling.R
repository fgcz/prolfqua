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
  res <- function(x){
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


#' get mixed model no interactions with peptide and sample random factor, from config
#' @export
model_no_interaction_and_sample_lmer <- function(config, factor_level = 2){
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
  broom::tidy(anova(modelNO,model))[2,"p.value"]
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


# extracting results ----

#' get all comparisons
#' @export
contrast_tukey_multcomp <- function(model, factor){
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
plot_lmer4_peptide_noRandom <- function(m){
  data <- m@frame
  ran <- ranef(m)[[1]]
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
  gg <- gg + theme(axis.text.x=element_text(angle = -90, hjust = 0))
  gg <- gg + geom_boxplot(alpha=0.1)
  return(gg)
}


#' get matrix of indicator coefficients for each interaction
#' @export
#'
lmer4_coeff_matrix <- function(m){
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
#' @export
coeff_weights_factor_levels <- function(mm){
  getCoeffs <- function(factor_level, mm){
    idx <- grep(factor_level, rownames(cm$mm))
    x<- as.list(apply(cm$mm[idx,],2,mean) )
    x <- as.tibble(x)
    add_column(x, "factor_level" = factor_level,.before=1)
  }
  factor_levels <- unique(unlist(str_split(rownames(mm), ":")))
  xx <- map_df(factor_levels, getCoeffs, mm)
  return(xx)
}




#' Add predicted values for each interaction
#' @export
#'
plot_predicted_interactions <- function(gg, m){
  cm <- lmer4_coeff_matrix(m)
  xstart_end <- data.frame(xstart = rownames(cm$mm), xend = rownames(cm$mm))
  ystart_end <- data.frame(xend = rownames(cm$mm), ystart =rep(0, nrow(cm$mm)),
                           yend = cm$mm %*% cm$coeffs)
  segments <- inner_join(xstart_end, ystart_end)
  gg <- gg + geom_segment(aes(x = xstart, y = ystart , xend = xend, yend =yend), data=segments, color = "blue", arrow=arrow())
  return(gg)
}

#' Make model plot with title - protein Name.
#' @export
#'
plot_model_and_data <- function(m, proteinID){
  gg <- plot_lmer4_peptide_noRandom(m)
  gg <- plot_predicted_interactions(gg, m)
  gg <- gg + ggtitle(proteinID)
  gg
}



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
#'
#' @export
workflow_lme4_model_analyse <- function(nestProtein, modelFunction, modelName, prot_stats)
{
  lmermodel <- paste0("lmer_", modelName)
  exists_lmer <- paste0("exists_lmer_", modelName)
  Coeffs_model <- paste0("Coeffs_", modelName)
  Anova_model <- paste0("Anova_", modelName)

  nestProtein %>% mutate(!!lmermodel := purrr::map(data, modelFunction)) ->
    modelProtein

  modelProtein <- modelProtein %>% mutate(!!exists_lmer := map_lgl(!!sym(lmermodel), function(x){!is.null(x)}))
  modelProteinF <- modelProtein %>% filter( !!sym(exists_lmer) == TRUE)
  no_ModelProtein <- modelProtein %>% filter(!!sym(exists_lmer) == FALSE)

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

  reslist$fig$VolcanoPlot <- Model_Coeff %>% filter(row.names.x. != "(Intercept)") %>%  quantable::multigroupVolcano(
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



