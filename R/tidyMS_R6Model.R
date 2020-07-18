#' build Model
#' @param data data - a data frame
#' @param modelFunction model function
#' @param grouping variable
#' @param modelName model name
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
#' mod <- LFQService:::build_model(
#'  pepIntensity,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = config$table$hkeysDepth())

build_model <- function(data,
                        modelFunction,
                        subject_Id = "protein_Id",
                        modelName = modelFunction$modelName){

  modellingResult <- LFQService:::model_analyse(
    pepIntensity,
    modelFunction,
    modelName = modelName,
    subject_Id = subject_Id)
  return( Model$new(modelDF = modellingResult$modelProtein,
                    modelName = modellingResult$modelName,
                    subject_Id = subject_Id))
}


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
#' mod <- LFQService:::build_model(
#'  pepIntensity,
#'  formula_randomPeptide,
#'  modelName = modelName,
#'  subject_Id = config$table$hkeysDepth())
#'
#' mod$modelDF
#' aovtable  <- mod$get_anova()
#' unique(aovtable$factor)
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



#' Do contrast
#' @export
#' @keywords internal
#' @examples
#'
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#'
#' D <- LFQServiceData::ionstar$normalized()
#' D$data <- dplyr::filter(D$data ,protein_Id %in% sample(protein_Id, 100))
#' modelName <- "f_condtion_r_peptide"
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
#'  mod$get_coefficients()
#'  Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")
#' contrast <- Contrasts$new(mod$modelDF, Contr, modelFunction, subject_Id = config$table$hkeysDepth())
#' contrast$get_contrasts_sides()
#' contrast$get_linfct()
#' contrast$get_contrasts()
#'
#' head(contrasts)
#' contrast$moderate()
#'
#' moderated_p_limma_long(contrasts, group_by_col = "contrast")
#'
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
    initialize = function(models, contrasts, modelFunction, subject_Id){
      self$models = models
      self$contrasts = contrasts
      self$contrastfun = modelFunction$contrast_fun
      self$modelName = modelFunction$model_name
      self$subject_Id = subject_Id
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
    #' @description get linear functions from contrasts
    get_linfct = function(){
      models <- self$models %>% dplyr::filter(exists_lmer == TRUE)
      m <- get_complete_model_fit(models)
      linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
      linfct <- unique(linfct) # needed for single factor models
      linfct_A <- linfct_matrix_contrasts(linfct, self$contrasts)
      return(list(linfct = linfct, linfct_A = linfct_A))
    },

    get_contrasts = function(){
      linfct <- self$get_linfct()
      contrast_sides <- self$get_contrasts_sides()

      contrast_result <- contrasts_linfct(self$models,
                                          rbind(linfct$linfct, linfct$linfct_A),
                                          subject_Id = self$subject_Id,
                                          contrastfun = self$contrastfun )
      contrast_result <- dplyr::rename(contrast_result, contrast = lhs)

      xx <- contrast_result %>% dplyr::select(self$subject_Id, "contrast", "estimate")
      xx <- xx %>% pivot_wider(names_from = "contrast", values_from = "estimate")

      contrast_result <- contrast_result %>% dplyr::filter(contrast %in% names(self$contrasts))

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
      return(contrast_result)
    },
    moderate = function(){
      res <- moderated_p_limma_long(self$get_contrasts(),group_by_col = "contrast")
      return(res)
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

# Contrasts ----
#' plot contrasts
#' @export
#' @examples
#'
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#'
#' D <- LFQServiceData::ionstar$normalized()
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
#'  mod$get_coefficients()
#'  Contr <- c("dil.b_vs_a" = "dilution.b - dilution.a",
#'   "dil.e_vs_a" = "dilution.e - dilution.a",
#'   "dil.e_vs_b" = "dilution.e - dilution.b",
#'   "dil.c_vs_b" = "dilution.c - dilution.b"
#'  )
#' contrast <- Contrasts$new(mod$modelDF, Contr, modelFunction, subject_Id = config$table$hkeysDepth())
#' #Contrast_Plotter$debug("volcano_plotly")
#' cp <- Contrast_Plotter$new(contrast$get_contrasts(),contrast$modelName, contrast$subject_Id)
#' cp$histogram()
#' cp$histogram_estimate()
#' res <- cp$volcano()
#' res[[1]]
#' res[[2]]
#' respltly <- cp$volcano_plotly()
#' cp$ma_plot()
#' cp$ma_plotly()
#' cp$write_pdf("c:/Temp")
#' cp$write_plotly("c:/Temp")
#' cp$get_minimal()
#' cp$to_wide()
#'
Contrast_Plotter <- R6::R6Class(
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
    #' @description create Crontrast_Plotter
    initialize = function(contrastDF,
                          modelName,
                          subject_Id,
                          estimate = "estimate",
                          contrast = "contrast"
                          ){
      self$contrastDF <- contrastDF %>% tidyr::unite("subject_Id", subject_Id, sep = "~", remove = FALSE)
      self$modelName  = modelName
      self$subject_Id = subject_Id
      self$estimate = estimate
    },
    #' @description  plot histogram of selected socres
    #' @param score which scores to show - list of lists
    histogram = function(score = list(list(score = "p.value", xlim = c(0,1,0.05)),
                                        list(score = "p.value.adjusted", xlim = c(0,1,0.05)),
                                        list(score = "statistic" , xlim = c(0,4,0.1)))
                         ){
      fig <- private$.histogram(score = score)

      self$figures[["histogram"]] <- list(fig = fig,
                                         name = paste0(self$prefix,"_Histogram_", self$modelName))
      return(fig)

    },
    #' @description plot histogram of estimated fold change
    histogram_estimate = function(){
      re <- range(self$contrastDF[[self$estimate]])
      re[1] <- floor(re[1])
      re[2] <- ceiling(re[2])
      fig <- private$.histogram(score = list(list(score =  self$estimate, xlim = c(re,0.1))))
      self$figures[["histogram_estimate"]] <- list(fig = fig,
                                         name = paste0(self$prefix,"_Histogram_Estimate_", self$modelName))
      return(fig)

    },
    #' @description plotly volcano plots
    #' @param scores for which scores to generate volcano plot
    #' @param  fc fold change abline
    volcano = function(scores = c("p.value","p.value.adjusted"), fc = 1){
      fig <- private$.volcano(self$contrastDF, scores, fc = fc )
      self$figures[["volcano"]] <- list(fig = fig, name = paste0(self$prefix, "_Volcano_", self$modelName))
      return(fig)
    },
    #' @description plotly volcano plots
    #' @param scores for which scores to generate volcano plot
    #' @param  fc fold change abline
    volcano_plotly = function(scores = c("p.value","p.value.adjusted"), fc = 1){
      contrastDF <- self$contrastDF %>% plotly::highlight_key(~ subject_Id)
      res <- private$.volcano(contrastDF, scores, fc = fc )
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

        self$figures_plotly[["ma_plot"]] <- list(fig = fig_plotly,
                                                name = paste0(self$prefix, "_MA_Plolty_", self$modelName))
        return(fig_plotly)

      }else{
        return(NULL)
      }
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
    write_pltly = function(path){
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
    },
    get_minimal = function(columns = c("p.value","p.value.adjusted")){
      relevant_columns <- c("contrast",
                            "c1_name",
                            "c1",
                            "c2_name",
                            "c2",
                            "sigma",
                            "df",
                            "isSingular",
                            "estimate",
                            "conf.low",
                            "conf.high") # other relevant columns.

      contrast_minimal <- self$contrastDF %>%
        dplyr::select_at(c(subject_Id, relevant_columns, columns) )
      return(contrast_minimal)

    },
    to_wide = function(columns = c("p.value","p.value.adjusted")){
      contrast_minimal <- self$get_minimal(columns = columns)
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = subject_Id,
                                                     columns = c("estimate", columns))
      return(contrasts_wide)
    }
  ),
  private = list(
    .volcano = function(contrasts, scores , fc = 1){
      fig <- list()
      for (column in scores) {
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
    .histogram  = function(score){
      plots <- list()
      for (i in score) {
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






