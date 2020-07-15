# Model -----

#' Model
#' @export
#' @examples
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#' D <- LFQServiceData::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' pepIntensity <- D$data
#' config <- D$config
#' config$table$hkeysDepth()
#' modellingResult <- LFQService:::model_analyse(
#'  pepIntensity,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = config$table$hkeysDepth())
#' names(modellingResult)
#'
#' mod <- Model$new(
#'   modelDF = modellingResult$modelProtein,
#'   modelName = modellingResult$modelName,
#'   subject_Id = config$table$hkeysDepth())
#'
#' mod$modelDF
#' mod$get_anova()
#' mod$get_coefficients()
#' mod$coef_histogram()
#' mod$coef_volcano()
#' mod$coef_pairs()
#' mod$anova_histogram()
Model <- R6::R6Class(
  "Model",
  public = list(
   #' @field modelDF data.frame with modelling data and model.
   modelDF = NULL,
   modelName = character(),
   subject_Id = character(),
   initialize = function(modelDF, modelName, subject_Id = "protein_Id"){
     self$modelDF = modelDF
     self$modelName = modelName
     self$subject_Id = subject_Id

   },
   get_coefficients = function(){
     lmermodel <- "linear_model"
     modelProteinF <- get_complete_model_fit(self$modelDF)
     # modelProteinF <- modelProteinF %>% dplyr::filter(nrcoef == max(nrcoef))
     # Extract coefficients
     .coef_df <-  function(x){
       x <- coef(summary(x));
       x <- data.frame(row.names(x), x);
       return(x)
     }
     Model_Coeff <- modelProteinF %>%
       dplyr::mutate(!!"Coeffs_model" := purrr::map( !!sym(lmermodel),  .coef_df ))
     Model_Coeff <- Model_Coeff %>%
       dplyr::select(!!!syms(self$subject_Id), !!sym("Coeffs_model"), isSingular, nrcoef)
     Model_Coeff <- tidyr::unnest_legacy(Model_Coeff)
     return(Model_Coeff)
   },
   get_anova = function(){
     lmermodel <- "linear_model"
     modelProteinF <- get_complete_model_fit(self$modelDF)
     # ANOVA
     .anova_df <- function(x){
       x <- anova(x)
       colnames(x) <- make.names(colnames(x))
       x <- data.frame(rownames(x), x)
       return(x)
     }

     Model_Anova <- modelProteinF %>% dplyr::mutate(!!"Anova_model" := purrr::map( !!sym(lmermodel),  .anova_df ))

     Model_Anova <- Model_Anova %>%
       dplyr::select(!!!syms(self$subject_Id), !!sym("Anova_model"), isSingular, nrcoef)
     Model_Anova <- tidyr::unnest_legacy(Model_Anova)
     #tidyr::unnest(cols = "Anova_model")
     return(Model_Anova)

   },


   #' writes results of `model_analyse`, anova table and all the coefficients with parameters.
   #' @keywords internal
   #' @export
   write_coefficients  = function(path){
       lfq_write_table(self$get_coefficients(),
                       path = path,
                       name  = paste0("Coef_",self$modelName) )
   },
   write_anova = function(path){
     lfq_write_table(self$get_anova(),
                     path = path,
                     name  = paste0("ANOVA_",self$modelName) )

   },

   #' @description histogram of model coefficient
   coef_histogram = function(){
     Model_Coeff <- self$get_coefficients()
     Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
     ## Coef_Histogram
     fname_histogram_coeff_p.values <- paste0("Coef_Histogram_",self$modelName,".pdf")
     histogram_coeff_p.values <- ggplot(data = Model_Coeff, aes(x = Pr...t.., group = row.names.x.)) +
       geom_histogram(breaks = seq(0,1,by = 0.05)) +
       facet_wrap(~row.names.x.)
     return(list(plot = histogram_coeff_p.values, name = fname_histogram_coeff_p.values))
   },
   #' @description volcano plot of non intercept coefficients
   coef_volcano = function(){
     Model_Coeff <- self$get_coefficients()
     Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
     fname_VolcanoPlot <- paste0("Coef_VolcanoPlot_",self$modelName,".pdf")
     VolcanoPlot <- Model_Coeff %>%
       dplyr::filter(row.names.x. != "(Intercept)") %>%
       LFQService::multigroupVolcano(
         effect = "Estimate",
         p.value = "Pr...t..",
         condition = "row.names.x.",
         label = "subject_Id" ,
         xintercept = c(-1, 1) ,
         colour = "isSingular" )
     return(list(plot = VolcanoPlot, name = fname_VolcanoPlot))
   },
   #' @description pairsplot of coefficients
   coef_pairs = function(){
     Model_Coeff <- self$get_coefficients()
     Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
     ## Coef_Pairsplot
     forPairs <- Model_Coeff %>%
       dplyr::select(!!sym("subject_Id") , row.names.x. ,  Estimate ) %>%
       tidyr::spread(row.names.x., Estimate )
     fname_Pairsplot_Coef <- paste0("Coef_Pairsplot_", self$modelName,".pdf")
     Pairsplot_Coef <-  GGally::ggpairs(forPairs, columns = 2:ncol(forPairs))
     return(list(plot = Pairsplot_Coef, name = fname_Pairsplot_Coef))

   },
   #' @description histogram of anova results
   anova_histogram = function(){
     ## Anova_p.values
     Model_Anova <- self$get_anova()
     fname_histogram_anova_p.values <- paste0("Anova_p.values_", self$modelName, ".pdf")
     histogram_anova_p.values <-  Model_Anova %>%
       dplyr::filter(rownames.x. != "Residuals") %>%
       ggplot( aes(x = Pr..F., group = rownames.x.)) +
       geom_histogram(breaks = seq(0,1,by = 0.05)) +
       facet_wrap(~rownames.x.)
     return(list(plot = histogram_anova_p.values, name = fname_histogram_anova_p.values))
   },
   write_anova_figures = function(path, width = 10, height =10){
     private$write_fig(self$anova_histogram(),path, width, height )
   },
   write_coef_figures = function(path, width = 10, height =10){
     private$write_fig(self$coef_histogram(),path, width, height )
     private$write_fig(self$coef_volcano(),path, width, height )
     private$write_fig(self$coef_pairs(),path, width, height )
   }
  ),
  private = list(
    write_fig = function(res, path, width = 10, height = 10){
      fpath <- file.path(path, res$name)
      message("Writing figure into : ", fpath, "\n")
      pdf(fpath, width = width, height = height )
      print(res$plot)
      dev.off()
    }
  )
)


#' p2621 workflow likelihood ratio test
#' @keywords internal
#' @export
#' @examples
#' #todo add example
workflow_likelihood_ratio_test <- function(modelProteinF,
                                           modelName,
                                           modelProteinF_Int,
                                           modelName_Int,
                                           subject_Id = "protein_Id",
                                           path = NULL
){
  # Model Comparison
  reg <- dplyr::inner_join(dplyr::select(modelProteinF, !!sym(subject_Id), "linear_model"),
                           dplyr::select(modelProteinF_Int, !!sym(subject_Id), "linear_model") , by = subject_Id)

  reg <- reg %>% dplyr::mutate(modelComparisonLikelihoodRatioTest = map2(!!sym("linear_model.x"),
                                                                         !!sym("linear_model.y"),
                                                                         .likelihood_ratio_test ))
  likelihood_ratio_test_result <- reg %>%
    dplyr::select(!!sym(subject_Id), modelComparisonLikelihoodRatioTest) %>%
    tidyr::unnest(cols = c("modelComparisonLikelihoodRatioTest"))
  likelihood_ratio_test_result <- likelihood_ratio_test_result %>%
    dplyr::rename(likelihood_ratio_test.pValue = modelComparisonLikelihoodRatioTest)


  if (!is.null(path)) {
    fileName <- paste("hist_LRT_", modelName, "_", modelName_Int, ".pdf", sep = "")
    fileName <- file.path(path, fileName)
    message("writing figure : " , fileName , "\n")
    pdf(fileName)
    par(mfrow = c(2, 1))
    hist(likelihood_ratio_test_result$likelihood_ratio_test.pValue,
         breaks = 20)
    plot(ecdf(
      likelihood_ratio_test_result$likelihood_ratio_test.pValue
    ))
    abline(v = c(0.01, 0.05), col = c(3, 2))
    dev.off()
  }

  return(likelihood_ratio_test_result)
}



