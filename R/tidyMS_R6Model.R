#' R6 class representing modelling result
#'
#' @export
#' @family modelling
#' @examples
#'
#' library(prolfqua)
#' library(tidyverse)
#' istar <- prolfqua::ionstar$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' pepIntensity <- istar_data
#' config <- istar$config
#' config$table$hkeysDepth()
#' mod <- prolfqua::build_model(
#'  pepIntensity,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = config$table$hkeysDepth())
#'
#' mod$modelDF
#' aovtable  <- mod$get_anova()
#' head(aovtable)
#' unique(aovtable$factor)
#' mod$get_coefficients()
#' mod$coef_histogram()
#' mod$coef_volcano()
#' mod$coef_pairs()
#' mod$anova_histogram()
#'
Model <- R6::R6Class(
  "Model",
  public = list(
    #' @field modelDF data.frame with modelling data and model.
    modelDF = NULL,
    #' @field modelName name of model
    modelName = character(),
    #' @field subject_Id e.g. protein_Id
    subject_Id = character(),
    #' @field modelFunction function to create the models
    modelFunction = NULL,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description initialize
    #' @param modelDF dataframe with modelling results
    #' @param modelFunction modelFunction see \code{\link{make_custom_model_lmer}}
    #' @param modelName name of model
    #' @param subject_Id subject column name
    #' @param p.adjust method to adjust p-values
    #'
    initialize = function(modelDF,
                          modelFunction,
                          modelName,
                          subject_Id = "protein_Id",
                          p.adjust = prolfqua::adjust_p_values){
      self$modelDF = modelDF
      self$modelFunction = modelFunction
      self$modelName = modelName
      self$subject_Id = subject_Id
      self$p.adjust = p.adjust
    },
    #' @description return model coefficient table
    get_coefficients = function(){
      lmermodel <- "linear_model"
      modelProteinF <- get_complete_model_fit(self$modelDF)
      # modelProteinF <- modelProteinF %>% dplyr::filter(nrcoef == max(nrcoef))
      # Extract coefficients
      .coef_df <-  function(x){
        x <- coef(summary(x));
        x <- data.frame(factor = row.names(x), x);
        return(x)
      }
      Model_Coeff <- modelProteinF %>%
        dplyr::mutate(!!"Coeffs_model" := purrr::map( !!sym(lmermodel),  .coef_df ))
      Model_Coeff <- Model_Coeff %>%
        dplyr::select(!!!syms(self$subject_Id), !!sym("Coeffs_model"), isSingular, nrcoef)
      Model_Coeff <- tidyr::unnest_legacy(Model_Coeff)
      return(Model_Coeff)
    },
    #' @description  return anova table
    get_anova = function(){
      lmermodel <- "linear_model"
      modelProteinF <- get_complete_model_fit(self$modelDF)
      # ANOVA
      .anova_df <- function(x){
        x <- anova(x, test = "F")
        colnames(x) <- make.names(colnames(x))
        x <- data.frame(factor = rownames(x), x)
        return(x)
      }

      Model_Anova <- modelProteinF %>% dplyr::mutate(!!"Anova_model" := purrr::map( !!sym(lmermodel),  .anova_df ))

      Model_Anova <- Model_Anova %>%
        dplyr::select(!!!syms(self$subject_Id), !!sym("Anova_model"), isSingular, nrcoef)
      Model_Anova <- tidyr::unnest_legacy(Model_Anova)



      Model_Anova <- Model_Anova %>% dplyr::filter(factor != "Residuals")
      Model_Anova <- Model_Anova %>% dplyr::filter(factor != "NULL")

      Model_Anova <- self$p.adjust(Model_Anova,
                                   column = "Pr..F.",
                                   group_by_col = "factor",
                                   newname = "FDR.Pr..F.")
      #tidyr::unnest(cols = "Anova_model")
      return(dplyr::ungroup(Model_Anova))
    },

    #' @description writes model coefficients to file
    #' @param path folder to write to
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write_coefficients  = function(path, format = "xlsx"){
      lfq_write_table(self$get_coefficients(),
                      path = path,
                      name  = paste0("Coef_",self$modelName),
                      format = format)
    },
    #' @description writes anova coefficients to file
    #' @param path folder to write to
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write_anova = function(path, format = "xlsx"){
      lfq_write_table(self$get_anova(),
                      path = path,
                      name  = paste0("ANOVA_",self$modelName) ,
                      format = format)

    },

    #' @description histogram of model coefficient
    coef_histogram = function(){
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      ## Coef_Histogram
      fname_histogram_coeff_p.values <- paste0("Coef_Histogram_",self$modelName,".pdf")
      histogram_coeff_p.values <- ggplot(data = Model_Coeff, aes(x = Pr...t.., group = factor)) +
        geom_histogram(breaks = seq(0,1,by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = histogram_coeff_p.values, name = fname_histogram_coeff_p.values))
    },
    #' @description volcano plot of non intercept coefficients
    coef_volcano = function(){
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      fname_VolcanoPlot <- paste0("Coef_VolcanoPlot_",self$modelName,".pdf")
      VolcanoPlot <- Model_Coeff %>%
        dplyr::filter(factor != "(Intercept)") %>%
        prolfqua::multigroupVolcano(
          effect = "Estimate",
          p.value = "Pr...t..",
          condition = "factor",
          label = "subject_Id" ,
          xintercept = c(-1, 1) ,
          colour = "isSingular" )
      return(list(plot = VolcanoPlot, name = fname_VolcanoPlot))
    },
    #' @description pairs-plot of coefficients
    coef_pairs = function(){
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      ## Coef_Pairsplot
      forPairs <- Model_Coeff %>%
        dplyr::select(!!sym("subject_Id") , factor ,  Estimate ) %>%
        tidyr::spread(factor, Estimate )
      fname_Pairsplot_Coef <- paste0("Coef_Pairsplot_", self$modelName,".pdf")
      Pairsplot_Coef <-  GGally::ggpairs(forPairs, columns = 2:ncol(forPairs))
      return(list(plot = Pairsplot_Coef, name = fname_Pairsplot_Coef))

    },
    #' @description histogram of ANOVA results
    #' @param what show either "Pr..F." or "FDR.Pr..F."
    anova_histogram = function(what=c("Pr..F.", "FDR.Pr..F.")){
      ## Anova_p.values
      what <- match.arg(what)
      Model_Anova <- self$get_anova()
      fname_histogram_anova_p.values <- paste0("Anova_p.values_", self$modelName, ".pdf")
      histogram_anova_p.values <-  Model_Anova %>%
        ggplot( aes(x = !!sym(what), group = factor)) +
        geom_histogram(breaks = seq(0,1,by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = histogram_anova_p.values, name = fname_histogram_anova_p.values))
    },
    #' @description write figures related to ANOVA into pdf file
    #' @param path folder name
    #' @param width figure width
    #' @param height figure height
    #'
    write_anova_figures = function(path, width = 10, height =10){
      private$write_fig(self$anova_histogram(),path, width, height )
    },
    #' @description write figures related to Coefficients into pdf file
    #' @param path folder name
    #' @param width figure width
    #' @param height figure height
    #'
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


#' likelihood ratio test
#' @keywords internal
#' @family modelling
#' @export
#' @examples
#' #todo add example
LR_test <- function(modelProteinF,
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


#' build from data and modelFunction Model
#'
#'
#'
#' @param data data - a data frame
#' @param modelFunction model function
#' @param subject_Id grouping variable
#' @param modelName model name
#' @return
#' a object of class \code{\link{Model}}
#' @family modelling
#' @seealso \code{\link{model_analyse}}, \code{\link{make_custom_model_lmer}} \code{\link{make_custom_model_lm}}
#'
#' @export
#' @examples
#' rm(list = ls())
#' library(prolfqua)
#' library(tidyverse)
#' D <- prolfqua::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#'
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' pepIntensity <- D$data
#' config <- D$config
#'
#'
#'
#' mod <- prolfqua:::build_model(
#'  pepIntensity,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = config$table$hkeysDepth())
#'
#'
build_model <- function(data,
                        modelFunction,
                        subject_Id = "protein_Id",
                        modelName = modelFunction$modelName){

  modellingResult <- prolfqua:::model_analyse(
    data,
    modelFunction,
    modelName = modelName,
    subject_Id = subject_Id)
  return( Model$new(modelDF = modellingResult$modelProtein,
                    modelFunction = modelFunction,
                    modelName = modellingResult$modelName,
                    subject_Id = subject_Id))
}
