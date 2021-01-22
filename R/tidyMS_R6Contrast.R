
# Contrast simple ----

#' compute contrasts with data imputation (directly from data)
#'
#' if interaction average can not be computed infer it using the 10\%
#'  smallest interaction averages in the dataset. Based on these compute fold changes.
#'  Use median of peptide level fold changes as protein estimate.
#'
#' @family modelling
#' @export
#' @examples
#'
#' library(LFQService)
#' library(tidyverse)
#' bb <- LFQService::ionstar$normalized()
#' configur <- bb$config
#' data <- bb$data
#' lfqdata <- LFQData$new(data, configur)
#' Contrasts <- c("dilution.b-a" = "dilution.b - dilution.a", "dilution.c-e" = "dilution.c - dilution.e")
#' tmp <- ContrastsSimpleImpute$new(lfqdata, contrasts = Contrasts)
#' tmp$get_contrast_sides()
#'
#' bb <- tmp$get_contrasts()
#' #bb$estimate_mad == 0
#' pl <- tmp$get_Plotter()
#' pl$histogram()
#' pl$histogram_estimate()
#' pl$ma_plot()
#' pl$volcano()
#'
ContrastsSimpleImpute <- R6::R6Class(
  "ContrastSimple",
  public = list(
    #' @field subject_Id subject_id e.g. protein_ID column
    subject_Id = character(),
    #' @field contrasts array with contrasts (see example)
    contrasts = character(),
    #' @field modelName model name
    modelName = character(),
    #' @field contrast_result data frame with results of contrast computation
    contrast_result = NULL,
    #' @description initialize
    #' @param lfqdata LFQData
    #' @param contrasts array of contrasts (see example)
    #' @param modelName default "groupAverage"
    initialize = function(lfqdata,
                          contrasts,
                          modelName = "groupAverage"){
      self$subject_Id = lfqdata$config$table$hkeysDepth()
      self$contrasts = contrasts
      self$modelName = modelName
      self$contrast_result = get_imputed_contrasts(lfqdata$data, lfqdata$config, contrasts)
      self$contrast_result$isSingular <- TRUE
      self$contrast_result <- mutate(self$contrast_result,
                                     statistic = estimate_median / estimate_mad )
      self$contrast_result <- mutate(self$contrast_result,
                                     p.value = 2*pt(abs(statistic) , df = 2 , lower.tail = FALSE) )
    },
    #' @description get contrasts sides
    #'
    get_contrast_sides = function(){
      self$contrast_result %>% select(contrast,c1 = c1_name,c2 = c2_name) %>% distinct()
    },
    #' @description table with results of contrast computation
    get_contrasts = function(){
      self$contrast_result
    },
    #' @description write contrasts computation results to
    #' @param path folder to write to
    #' @param format default xlsx (can be csv or html)
    write = function(path, format = "xlsx"){
      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",self$modelName) ,
                      format = format)
    },
    #' @description get Contrast_Plotter
    get_Plotter = function(){
      res <- Contrasts_Plotter$new(
        self$contrast_result,
        subject_Id = self$subject_Id,
        volcano = list(list(score = "p.value", fc = 1)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "estimate_median",
        contrast = "contrast")
      return(res)
    },
    #' @description convert contrast results to wide format
    #' @param columns value column default p.value
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("estimate_median"),
                                                     contrast = 'contrast')
      return(contrasts_wide)
    }

  )
)



# Contrasts -----

#' Estimate contrasts using Wald Test
#' @export
#' @family modelling
#' @examples
#'
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#'
#' istar <- LFQService::ionstar$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id) + (1 | sampleName)")
#' pepIntensity <- istar$data
#' config <- istar$config
#' config$table$hkeysDepth()
#'
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  subject_Id = config$table$hkeysDepth())
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b",
#'     "dil.e_vs_a" = "dilution.e - dilution.b" )
#'  #LFQService::Contrasts$debug("get_linfct")
#' contrastX <- LFQService::Contrasts$new(mod,
#'  Contr)
#'
#' contrastX$get_contrasts_sides()
#'
#' contrastX$get_linfct()
#' xx <- contrastX$get_contrasts()
#'
#'
#'
Contrasts <- R6::R6Class(
  "Contrast",
  public = list(
    #' @field models Model
    models = NULL,
    #' @field contrasts character with contrasts
    contrasts = character(),
    #' @field contrastfun function to compute contrasts
    contrastfun = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id name of column containing e.g., protein Id's
    subject_Id = character(),
    #' @field p.adjust function to adjust p-values (default LFQService::adjust_p_values)
    p.adjust = NULL,
    #' @field contrast_result data frame containing results of contrast computation
    contrast_result = NULL,
    #' @description initialize
    #' create Contrast
    #' @param model a dataframe with a structure similar to that generated by \code{\link{build_model}}
    #' @param contrasts a character vector with contrast specificiation
    #' @param p.adjust function to adjust the p-values
    initialize = function(model,
                          contrasts,
                          p.adjust = LFQService::adjust_p_values
    ){
      self$models = model$modelDF
      self$contrasts = contrasts
      self$contrastfun = model$modelFunction$contrast_fun
      self$modelName =  model$modelFunction$model_name
      self$subject_Id = model$subject_Id
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
    #' @description get linear functions from contrasts
    #' @param global logical TRUE - get the a linear functions for all models, FALSE - linear function for each model
    get_linfct = function(global = TRUE){
      linfct <- function(model, contrast){
        linfct <- linfct_from_model(model, as_list = FALSE)
        linfct <- unique(linfct) # needed for single factor models
        linfct_A <- linfct_matrix_contrasts(linfct, self$contrasts)
        return( rbind( linfct, linfct_A ) )
      }
      if (global) {
        models <- self$models %>% dplyr::filter(exists_lmer == TRUE)
        model <- get_complete_model_fit(models)$linear_model[[1]]
        return( linfct( model, self$contrasts ))
      }else{
        linfct <- purrr::map(self$models$linear_model,
                             linfct, contrast = self$contrast)
        return(linfct)
      }
    },
    #' @description
    #' get table with contrast estimates
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    #' @return data.frame with contrasts
    get_contrasts = function(all = FALSE, global = TRUE){
      if (!is.null(self$contrast_result) ) {
        return(self$contrast_result)
      }

      message("determine linear functions:")
      linfct <- self$get_linfct(global = global)
      contrast_sides <- self$get_contrasts_sides()
      message("compute contrasts:")
      contrast_result <- contrasts_linfct(self$models,
                                          linfct,
                                          subject_Id = self$subject_Id,
                                          contrastfun = self$contrastfun )

      contrast_result <- dplyr::rename(contrast_result, contrast = lhs)

      xx <- dplyr::select(contrast_result, self$subject_Id, "contrast", "estimate")
      xx <- xx %>% tidyr::pivot_wider(names_from = "contrast", values_from = "estimate")

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
      contrast_result <- inner_join(contrast_sides, contrast_result)
      if (!all) {
        contrast_result <- contrast_result %>% select(-all_of(c("sigma.model",
                                                                "df.residual.model")))
      }
      self$contrast_result <- self$p.adjust(contrast_result,
                                            column = "p.value",
                                            group_by_col = "contrast",
                                            newname = "FDR")
      self$contrast_result <- dplyr::ungroup(self$contrast_result )
      return(self$contrast_result)
    },
    #' @description write results
    #' @param path directory
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write = function(path, format = "xlsx"){
      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",self$modelName),
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
    #' @description convert to wide format
    #' @param columns value column default p.value
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
#' @family modelling
#' @examples
#'
#' library(LFQService)
#' istar <- LFQService::ionstar$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id) + (1|sampleName)")
#' pepIntensity <- istar_data
#' config <- istar$config$clone(deep = TRUE)
#' config$table$hkeysDepth()
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  subject_Id = config$table$hkeysDepth())
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")
#'  contrast <- LFQService::ContrastsModerated$new(mod,
#'  Contr)
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
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE, global = TRUE){
      contrast_result <- super$get_contrasts(all = TRUE, global = global)
      contrast_result <- moderated_p_limma_long(contrast_result , group_by_col = "contrast")
      if (!all) {
        contrast_result <- select(contrast_result ,
                                  -all_of(c("sigma.model",
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
    #' @description get \code{\link{Contrast_Plotter}}
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
    #' @description convert to wide format
    #' @param columns value column default moderated.p.value
    to_wide = function(columns = c("moderated.p.value", "FDR.moderated")){
      super$to_wide(columns = columns)
    }
  )
)


# ContrastsROPECA -----

#'
#' ROPECA reproducibility-optimization method
#'
#' ROPECA optimizes the reproducibility of statistical testing
#'  on peptide-level and aggregates the peptide-level changes
#'   to determine differential protein-level expression.
#'
#' @export
#' @family modelling
#' @examples
#'
#' istar <- LFQService::ionstar$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   make_custom_model_lm("transformedIntensity  ~ dilution.")
#' pepIntensity <- istar_data
#' config <- istar$config$clone(deep = TRUE)
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
#'  contrast <- LFQService::ContrastsROPECA$new(mod, Contr)
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
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all=FALSE, global = TRUE){
      contrasts_data <- super$get_contrasts(all = all, global = global)
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
    #' @description get \code{\link{Contrast_Plotter}}
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
    #' @description convert to wide format
    #' @param columns value column default beta.based.significance
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

  fig <- ggpubr::annotate_figure(fig, bottom = ggpubr::text_grob(annot, size = 10))
  return(fig)
}

# Contrast_Plotter ----
#' plot contrasts
#' @export
#' @family modelling
#' @family plotting
#' @examples
#'
#' rm(list = ls())
#' library(LFQService)
#' library(tidyverse)
#'
#' istar <- LFQService::ionstar$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelName <- "Model"
#' modelFunction <-
#'   make_custom_model_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
#'    model_name = modelName)
#' pepIntensity <- istar_data
#' config <- istar$config
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
#' contrast <- LFQService::Contrasts$new(mod,
#'   Contr)
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
#' cp$ma_plotly
#' if(FALSE){
#'   cp$write_pdf("c:/Temp")
#'   cp$write_plotly("c:/Temp")
#' }
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
    #' @field estimate column with fold change estimates
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
    #' @param subject_Id columns containing subject Identifier
    #' @param volcano which score to plot and which ablines to add.
    #' @param histogram which scores to plot and which range (x) should be shown.
    #' @param modelName name of model
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
    #' @description  plot histogram of selected scores (e.g. p-value, FDR, t-statistics)
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
    #' @description plotly volcano plots (fold change vs FDR)
    volcano = function(){
      fig <- private$.volcano(self$contrastDF, self$volcano_spec )
      self$figures[["volcano"]] <- list(fig = fig, name = paste0(self$prefix, "_Volcano_", self$modelName))
      return(fig)
    },
    #' @description plotly volcano plots
    #' @param scores for which scores to generate volcano plot
    #' @param fc fold change abline
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
    #' @description write figures in pdf format
    #' @param path path
    #' @param width figure with
    #' @param height figure height
    #'
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
    #' @description write figures into thml
    #' @param path path
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
        fig[[column]] <- LFQService:::.multigroupVolcano(
          contrasts,
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






