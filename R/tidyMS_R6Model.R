

# Model -----

#' Do Model
#'
#' @export
#' @examples
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#' D <- LFQService::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#' modelName <- "f_condtion_r_peptide"
#' formula_randomPeptide <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' pepIntensity <- D$data
#' config <- D$config
#' config$table$hkeysDepth()
#' mod <- LFQService::build_model(
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
    modelName = character(),
    subject_Id = character(),
    initialize = function(modelDF, modelName = "Model", subject_Id = "protein_Id"){
      self$modelDF = modelDF
      self$modelName = modelName
      self$subject_Id = subject_Id

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
    #' return anova table
    get_anova = function(){
      lmermodel <- "linear_model"
      modelProteinF <- get_complete_model_fit(self$modelDF)
      # ANOVA
      .anova_df <- function(x){
        x <- anova(x)
        colnames(x) <- make.names(colnames(x))
        x <- data.frame(factor = rownames(x), x)
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
    write_coefficients  = function(path, format = "xlsx"){
      lfq_write_table(self$get_coefficients(),
                      path = path,
                      name  = paste0("Coef_",self$modelName),
                      format = format)
    },
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
        LFQService::multigroupVolcano(
          effect = "Estimate",
          p.value = "Pr...t..",
          condition = "factor",
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
        dplyr::select(!!sym("subject_Id") , factor ,  Estimate ) %>%
        tidyr::spread(factor, Estimate )
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
        dplyr::filter(factor != "Residuals") %>%
        ggplot( aes(x = Pr..F., group = factor)) +
        geom_histogram(breaks = seq(0,1,by = 0.05)) +
        facet_wrap(~factor)
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
#' @family modelling
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


#' build from data and modelFunction Model
#'
#'
#'
#' @param data data - a data frame
#' @param modelFunction model function
#' @param grouping variable
#' @param modelName model name
#' @return
#' a object of class \code{\link{Model}}
#' @family modelling
#' @seealso \code{\link{model_analyse}}, \code{\link{make_custom_model_lmer}} \code{\link{make_custom_model_lm}}
#'
#' @export
#' @examples
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#' D <- LFQService::ionstar$normalized()
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
#' mod <- LFQService:::build_model(
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

  modellingResult <- LFQService:::model_analyse(
    data,
    modelFunction,
    modelName = modelName,
    subject_Id = subject_Id)
  return( Model$new(modelDF = modellingResult$modelProtein,
                    modelName = modellingResult$modelName,
                    subject_Id = subject_Id))
}


# Contrasts -----

#' Do contrast
#' @export
#' @examples
#'
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#'
#' D <- LFQService::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id) + (1 | sampleName)")
#' pepIntensity <- D$data
#' config <- D$config
#' config$table$hkeysDepth()
#'
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  subject_Id = config$table$hkeysDepth())
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")
#'  #LFQService::Contrasts$debug("get_linfct")
#' contrastX <- LFQService::Contrasts$new(mod$modelDF,
#'  Contr,
#'  modelFunction$contrast_fun,
#'  subject_Id = config$table$hkeysDepth(),
#'  modelName = modelFunction$model_name)
#'
#' contrastX$get_contrasts_sides()
#'
#' contrastX$get_linfct()
#' xx <- contrastX$get_contrasts()
#'
#'
#' imputed <- get_imputed_contrasts(D$data, D$config, Contr)
#'
Contrasts <- R6::R6Class(
  "Contrast",
  public = list(
    #' @field models Model
    models = NULL,
    #' @field contrasts character with contrasts
    contrasts = character(),
    #' @field contrast function
    contrastfun = NULL,
    #' @field modelName model name
    modelName = NULL,
    subject_Id = character(),
    #' @field function to adjust p-values (default LFQService::adjust_p_values)
    p.adjust = NULL,
    contrast_result = NULL,
    #' create Contrast
    #' @param models a dataframe with a structure similar to that generated by \code{\link{build_model}}
    #' @param contrasts a character vector with contrast specificiation
    #' @param contrast_fun see. eg. \code{\link{my_contrast_V2},\link{my_contest}}
    #' @param subject_Id columns with subject_Id (e.g. protein_Id)
    #' @param modelName default "Model"
    #' @param p.adjust function to adjust the p-values
    #' @param contrast_result contrast result
    initialize = function(models,
                          contrasts,
                          contrast_fun,
                          subject_Id,
                          modelName = "Model",
                          p.adjust = LFQService::adjust_p_values
                          ){
      self$models = models
      self$contrasts = contrasts
      self$contrastfun = contrast_fun
      self$modelName = modelName
      self$subject_Id = subject_Id
      self$p.adjust = p.adjust
    },
    #' @description get both sides of contrasts
    get_contrasts_sides = function(){
      # extract contrast sides
      tt <- self$contrasts[grep("-",self$contrasts)]
      tt <- tibble(contrast = names(tt) , rhs = tt)
      tt <- tt %>% mutate(rhs = gsub("[` ]","",rhs)) %>%
        tidyr::separate(rhs, c("c1", "c2"), sep = "-")
      return(tt)
    },
    #' @description
    #' get linear functions from contrasts
    get_linfct = function(){
      models <- self$models %>% dplyr::filter(exists_lmer == TRUE)
      m <- get_complete_model_fit(models)
      linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
      linfct <- unique(linfct) # needed for single factor models
      linfct_A <- linfct_matrix_contrasts(linfct, self$contrasts)
      return(list(linfct = linfct, linfct_A = linfct_A))
    },
    #' @description
    #' get table with contrast estimates
    #' @param all should all columns be returned (default FALSE)
    #' @return data.frame with contrasts
    get_contrasts = function(all = FALSE){
      if (!is.null(self$contrast_result) ) {
        return(self$contrast_result)
      }
      linfct <- self$get_linfct()
      contrast_sides <- self$get_contrasts_sides()
      contrast_result <- contrasts_linfct(self$models,
                                          rbind(linfct$linfct, linfct$linfct_A),
                                          subject_Id = self$subject_Id,
                                          contrastfun = self$contrastfun )
      contrast_result <- dplyr::rename(contrast_result, contrast = lhs)

      xx <- dplyr::select(contrast_result, self$subject_Id, "contrast", "estimate")
      xx <- xx %>% pivot_wider(names_from = "contrast", values_from = "estimate")

      contrast_result <- dplyr::filter(contrast_result, contrast %in% names(self$contrasts))

      get_contrast_cols <- function(i, contrast_results , contrast_table , subject_ID ){
        data.frame(lhs = contrast_table[i, "contrast"],
                   dplyr::select_at(contrast_results, c( subject_ID, unlist(contrast_table[i,c("c1","c2")]))),
                   c1_name = contrast_table[i,"c1", drop = T],
                   c2_name = contrast_table[i,"c2", drop = T], stringsAsFactors = FALSE)
      }

      contrast_sides <- purrr::map_df(1:nrow(contrast_sides),
                                      get_contrast_cols,
                                      xx,
                                      contrast_sides,
                                      self$subject_Id)
      contrast_result <- inner_join(contrast_sides,contrast_result)
      if (!all) {
        contrast_result <- contrast_result %>% select(-all_of(c("sigma.model",
                                                                "df.residual.model")))
      }

      self$contrast_result <- self$p.adjust(contrast_result,
                                  column = "p.value",
                                  group_by_col = "contrast",
                                  newname = "FDR")
      return(self$contrast_result)
    },
    write = function(path, format = "xlsx"){
      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",self$modelName) ,
                      format = format)
    },
    #' @description
    #' return Contrast_Plotter
    get_Plotter = function(){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        volcano = list(list(score = "FDR", fc = 1)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "estimate",
        contrast = "contrast")
      return(res)
    },
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("estimate", columns),
                                                     contrast = 'contrast')
      return(contrasts_wide)
    }
  ))

# ContrastsModerated -----


#' Limma moderated contrasts
#' @export
#'
#' @examples
#' D <- LFQService::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id) + (1|sampleName)")
#' pepIntensity <- D$data
#' config <- D$config$clone(deep = TRUE)
#' config$table$hkeysDepth()
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  subject_Id = config$table$hkeysDepth())
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")
#'  contrast <- LFQService::ContrastsModerated$new(mod$modelDF,
#'  Contr,
#'  modelFunction$contrast_fun,
#'  subject_Id = config$table$hkeysDepth(),
#'  modelName = modelFunction$model_name)
#'  contrast$get_contrasts(all = TRUE)
#'  plotter <- contrast$get_Plotter()
#'  plotter$histogram()
#'  plotter$ma_plot()
#'  plotter$volcano()
#'
ContrastsModerated <- R6::R6Class(
  classname = "ContrastsModerated",
  inherit = Contrasts,
  public = list(
    #' @description applies limma moderation
    #' @seealso \code{\link{moderated_p_limma_long}}
    get_contrasts = function(all = FALSE){
      contrast_result <- super$get_contrasts(all = TRUE)
      contrast_result <- moderated_p_limma_long(contrast_result , group_by_col = "contrast")
      if (!all) {
        contrast_result <- select(contrast_result ,-all_of(c("sigma.model",
                                        "df.residual.model",
                                        "moderated.df.prior",
                                        "moderated.var.prior",
                                        "moderated.var.post",
                                        "moderated.statistic",
                                        "moderated.df.total" )))
      }
      contrast_result <- self$p.adjust(contrast_result,
                                       column = "moderated.p.value",
                                       group_by_col = "contrast",
                                       newname = "FDR.moderated")
      return(contrast_result)
    },
    get_Plotter = function(){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        volcano = list(list(score = "FDR.moderated", fc = 1)),
        histogram = list(list(score = "moderated.p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR.moderated", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "estimate",
        contrast = "contrast")
      return(res)
    },
    to_wide = function(columns = c("moderated.p.value", "FDR.moderated")){
      super$to_wide(columns = columns)
    }
  )
)


# ContrastsROPECA -----

#'
#' Ropeca
#'
#' @export
#' @examples
#' D <- LFQService::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   make_custom_model_lm("transformedIntensity  ~ dilution.")
#' pepIntensity <- D$data
#' config <- D$config$clone(deep = TRUE)
#' config$table$hierarchyDepth <- 2
#' config$table$hkeysDepth()
#'
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  subject_Id = config$table$hkeysDepth())
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")
#'
#'
#'  #ContrastsROPECA$debug("get_Plotter")
#'  contrast <- LFQService::ContrastsROPECA$new(mod$modelDF,
#'  Contr,
#'  modelFunction$contrast_fun,
#'  subject_Id = config$table$hkeysDepth(),
#'  modelName = modelFunction$model_name)
#'
#'  contrast$get_linfct()
#'  contrast$subject_Id
#'  tmp <- contrast$get_contrasts(all = TRUE)
#'
#'  pl <- contrast$get_Plotter()
#'
#'  pl$histogram()
#'  pl$ma_plot()
#'
ContrastsROPECA <- R6::R6Class(
  "ContrastsROPECA",
  inherit = ContrastsModerated,
  public = list(
    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    get_contrasts = function(all=FALSE){
      contrasts_data <- super$get_contrasts(all = TRUE)
      res <- summary_ROPECA_median_p.scaled(
        contrasts_data,
        contrast = "contrast",
        subject_Id = self$subject_Id[1],
        estimate = "estimate",
        statistic = "statistic",
        p.value = "moderated.p.value",
        max.n = 10)

      if (!all) {
        res <- select(res , -all_of(c( "n_not_na", "mad.estimate", "n.beta")) )
      }

      res <- self$p.adjust(res,
                           column = "beta.based.significance",
                           group_by_col = "contrast",
                           newname = "FDR.beta.based.significance")
      res <- self$p.adjust(res,
                           column = "median.p.value",
                           group_by_col = "contrast",
                           newname = "FDR.median.p.value")
      return(res)
    },
    get_Plotter = function(){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id[1],
        volcano = list(list(score = "FDR.beta.based.significance", fc = 1)),
        histogram = list(list(score = "beta.based.significance", xlim = c(0,1,0.05)),
                         list(score = "FDR.beta.based.significance", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "estimate",
        contrast = "contrast")
      return(res)
    },
    to_wide = function(columns = c("beta.based.significance", "FDR.beta.based.significance")){
      super$to_wide(columns = columns)
    }
  ))



# plot score distributions by species
.plot_score_distribution <- function(data,
                                     score =list(list(score = "estimate",xlim = c(-1,2) ),
                                                 list(score = "statistic", xlim = c(-3,10) )),
                                     contrast = "contrast",

                                     annot = "peptide level statistics density"){
  plots <- list()
  for (i in score) {
    xlim = i$xlim
    score = i$score
    plots[[score]] <- data %>% ggplot(aes(x = !!sym(score),
                                          y = !!sym(contrast))) +
      ggridges::geom_density_ridges(alpha = 0.1) + xlim(xlim)
  }
  fig <- ggpubr::ggarrange(plotlist = plots,
                           nrow = 1,
                           common.legend = TRUE,
                           legend = "bottom")

  fig <- ggpubr::annotate_figure(fig, bottom = text_grob(annot, size = 10))
  return(fig)
}

# Contrast_Plotter ----
#' plot contrasts
#' @export
#' @examples
#'
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#'
#' D <- LFQService::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#' modelName <- "Model"
#' modelFunction <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' pepIntensity <- D$data
#' config <- D$config
#' config$table$hkeysDepth()
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  modelName = modelName,
#'  subject_Id = config$table$hkeysDepth())
#'
#'  #mod$get_coefficients()
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.b - dilution.a",
#'   "dil.e_vs_a" = "dilution.e - dilution.a",
#'   "dil.e_vs_b" = "dilution.e - dilution.b",
#'   "dil.c_vs_b" = "dilution.c - dilution.b"
#'  )
#' #Contrasts$debug("get_contrasts")
#' contrast <- LFQService::Contrasts$new(mod$modelDF,
#'   Contr,
#'   modelFunction$contrast_fun,
#'   subject_Id = config$table$hkeysDepth(),
#'   modelName = modelFunction$modelName)
#' tmp <- contrast$get_contrasts()
#'
#' cp <- Contrasts_Plotter$new(tmp ,
#'  contrast$subject_Id,
#' volcano = list(list(score = "FDR", fc = 1)),
#' histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
#'                  list(score = "FDR", xlim = c(0,1,0.05))))
#' cp$histogram()
#' cp$histogram_estimate()
#' res <- cp$volcano()
#' length(res)
#' res
#' respltly <- cp$volcano_plotly()
#' length(respltly)
#' cp$ma_plot()
#' cp$ma_plotly()
#' cp$write_pdf("c:/Temp")
#' cp$write_plotly("c:/Temp")
#'
Contrasts_Plotter <- R6::R6Class(
  "Contrast_Plotter"
  ,
  public = list(
    #' @field contrastDF data frame with contrasts
    contrastDF = NULL,
    #' @field modelName name of model
    modelName =  character(),
    #' @field subject_Id hierarchy key columns
    subject_Id = character(),
    #' @field prefix default Contrasts - used to generate file names
    prefix = "Contrasts",
    #' column with fold change estimates
    estimate = "estimate",
    #' @field contrast column with contrasts names, default "contrast"
    contrast = "contrast",
    #' @field figures list of generated figures
    figures = list(),
    #' @field figures_plotly list of generated figures
    figures_plotly = list(),
    #' @field volcano_spec volcano plot specification
    #'
    volcano_spec = NULL,
    #' @field histogram_spec plot specification
    #'
    histogram_spec = NULL,
    #' @description create Crontrast_Plotter
    #' @param contrastDF frame with contrast data
    #' @param modelName name of model
    #' @param subject_Id columns containing subject Identifier
    #' @param volcano which score to plot and which ablines to add.
    #' @param histogram which scores to plot and which range (x) should be shown.
    #' @param modelName Model
    #' @param estimate estimate column
    #' @param contrast contrast column
    initialize = function(contrastDF,
                          subject_Id,
                          volcano = list(list(score = "FDR", fc = 1)),
                          histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                                           list(score = "FDR", xlim = c(0,1,0.05))),
                                           #list(score = "statistic" , xlim = c(0,4,0.1))),
                          modelName = "Model",
                          estimate = "estimate",
                          contrast = "contrast"
    ){
      self$contrastDF <- contrastDF %>% tidyr::unite("subject_Id", subject_Id, sep = "~", remove = FALSE)
      self$modelName  = modelName
      self$subject_Id = subject_Id
      self$estimate = estimate
      self$volcano_spec = volcano
      self$histogram_spec = histogram

    },
    #' @description  plot histogram of selected socres
    #' @param score which scores to show - list of lists
    histogram = function(){
      if (length(self$histogram_spec) > 0) {
        fig <- private$.histogram(scores = self$histogram_spec)

        self$figures[["histogram"]] <-
          list(fig = fig,
               name = paste0(self$prefix,"_Histogram_", self$modelName))
        return(fig)
      }
    },
    #' @description plot histogram of estimated fold change
    histogram_estimate = function(){
      re <- range(self$contrastDF[[self$estimate]], na.rm = TRUE)
      re[1] <- floor(re[1])
      re[2] <- ceiling(re[2])
      fig <- private$.histogram(score = list(list(score =  self$estimate, xlim = c(re,0.1))))
      self$figures[["histogram_estimate"]] <- list(fig = fig,
                                                   name = paste0(self$prefix,"_Histogram_Estimate_", self$modelName))
      return(fig)

    },
    #' @description plotly volcano plots
    volcano = function(){
      fig <- private$.volcano(self$contrastDF, self$volcano_spec )
      self$figures[["volcano"]] <- list(fig = fig, name = paste0(self$prefix, "_Volcano_", self$modelName))
      return(fig)
    },
    #' @description plotly volcano plots
    #' @param scores for which scores to generate volcano plot
    #' @param  fc fold change abline
    volcano_plotly = function(){
      contrastDF <- self$contrastDF %>% plotly::highlight_key(~ subject_Id)
      res <- private$.volcano(contrastDF, self$volcano_spec )
      for (i in 1:length(res)) {
        res[[i]] <- res[[i]] %>% plotly::ggplotly(tooltip = "subject_Id")
      }
      self$figures_plotly[["volcano"]] <- list(fig = res,
                                               name = paste0(self$prefix, "_Volcano_Plolty_", self$modelName))
      return(res)
    },
    #' @description
    #' ma plot
    #' @param fc fold change abline
    ma_plot = function(fc = 1){
      contrastDF <- self$contrastDF
      if (!is.null(contrastDF$c1) && !is.null(contrastDF$c2)) {
        # pdf version
        fig <- private$.ma_plot(contrastDF ,self$contrast, fc = fc)
        self$figures[["ma_plot"]] <- list(fig = fig,
                                          name = paste0(self$prefix, "_MA_", self$modelName))
      }else{
        warning("no c1 c2 columns can't generate MA")
        fig <- NULL
      }
      return(fig)
    },
    #' @description
    #' ma plotly
    #' @param fc fold change abline
    ma_plotly = function(fc =1){
      # html version
      contrastDF  <- self$contrastDF
      if (!is.null(contrastDF$c1) && !is.null(contrastDF$c2)) {

        contrastDF  <- contrastDF %>%
          plotly::highlight_key(~subject_Id)
        fig_plotly <- private$.ma_plot(contrastDF, self$contrast, fc = fc) %>%
          plotly::ggplotly(tooltip = "subject_Id")

        self$figures_plotly[["ma_plot"]] <-
          list(fig = fig_plotly,
               name = paste0(self$prefix, "_MA_Plolty_", self$modelName))
        return(fig_plotly)

      }else{
        return(NULL)
      }
    },
    #' @description
    #' generate all figures
    #'
    all_figs = function(){
      self$ma_plotly()
      self$volcano_plotly()
      self$ma_plot()
      self$volcano()
      self$histogram()
      self$histogram_estimate()
      invisible(NULL)
    },
    write_pdf = function(path,
                         width = 10,
                         height = 10){
      message("Writing: ",path,"\n")

      plotX <- function(fig, width, height, path, fname, xname){
        fpath <- file.path(path,
                           paste(c(fname, if (!is.null(xname)) { c("_" , xname) }, ".pdf"),
                                 collapse = ""))
        pdf(fpath, width = width, height = height)
        print(fig)
        dev.off()
      }

      for (fig in self$figures) {
        if ("list" %in% class(fig$fig)) {
          for (i in 1:length(fig$fig)) {
            plotX(fig$fig[[i]], width, height, path, fname = fig$name, xname = names(fig$fig)[i] )
          }
        }else{
          plotX(fig$fig, width, height, path, fname = fig$name, xname = NULL )
        }

      }
    },
    write_plotly = function(path){
      message("Writing: ",path,"\n")

      plotX <- function(fig,path, fname, xname){
        fname <- paste(c(fname, if (!is.null(xname)) { c("_" , xname) }, ".html"),
                       collapse = "")
        fpath <- file.path(path,fname)
        message("Writing: ",fpath,"\n")
        htmlwidgets::saveWidget(widget = fig, fname, selfcontained = TRUE)
        file.rename(fname, fpath)
      }

      for (fig in self$figures_plotly) {
        if ("list" %in% class(fig$fig)) {
          for (i in 1:length(fig$fig)) {
            plotX(fig$fig[[i]], path, fname = fig$name, xname = names(fig$fig)[i] )
          }
        }else{
          plotX(fig$fig,  path, fname = fig$name, xname = NULL )
        }

      }
    }
  ),
  private = list(
    .volcano = function(contrasts, scores){
      fig <- list()
      for (score in scores) {
        column <- score$score
        fc <- score$fc
        fig[[column]] <- LFQService:::.multigroupVolcano(contrasts,
                                                         effect = self$estimate,
                                                         p.value = column,
                                                         condition = self$contrast,
                                                         text = "subject_Id",
                                                         xintercept = c(-fc, fc),
                                                         colour = "isSingular",
                                                         scales = "free_y")

      }
      return(fig)

    },
    .histogram  = function(scores){
      plots <- list()
      for (i in scores) {
        xlim = i$xlim
        score = i$score
        plots[[score]] <- self$contrastDF %>% ggplot(aes(x = !!sym(score))) +
          geom_histogram(breaks = seq(from = xlim[1], to = xlim[2], by = xlim[3])) +
          facet_wrap(vars(!!sym(self$contrast)))
      }
      fig <- ggpubr::ggarrange(plotlist = plots,
                               nrow = 1,
                               common.legend = TRUE,
                               legend = "bottom")
      #annot <- "histogram of score distribution"
      #fig <- ggpubr::annotate_figure(fig, bottom = ggpubr::text_grob(annot, size = 10))
      return(fig)
    },
    .ma_plot = function(x, contrast,  fc = 1){
      x <- ggplot(x , aes(x = (c1 + c2)/2,
                          y = !!sym(self$estimate),
                          text = !!sym("subject_Id"),
                          colour = !!sym("isSingular"))) +
        geom_point(alpha = 0.5) +
        scale_colour_manual(values = c("black", "red")) +
        facet_wrap(vars(!!sym(self$contrast))) + theme_light() +
        geom_hline(yintercept = c(-fc, fc), linetype = "dashed",colour = "red")
      return(x)
    }


  )
)






