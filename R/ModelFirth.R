# Model -----
#' R6 class representing modelling result
#'
#' @export
#' @family modelling
#' @examples
#'
#'
#'
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
#' istar$data <- prolfqua::encode_bin_resp(istar$data, istar$config)
#' tmp <- as_lfq(istar)
#' formula <- paste0(tmp$config$table$bin_resp , "~ group_")
#' mod <- build_model_logistf(tmp, formula)
#' tmp <- mod$get_coefficients()
#' mod$coef_histogram()
#' mod$coef_pairs()
#' mod$get_anova()
#' mod$coef_volcano()
#' mod$anova_histogram()
#' mod$write_coef_figures(tempdir())
#'
#' istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 10, with_missing = TRUE, weight_missing = 0.5, seed = 3)
#' istar$data <- prolfqua::encode_bin_resp(istar$data, istar$config)
#' tmp <- as_lfq(istar)
#' formula <- paste0(tmp$config$table$bin_resp , "~ group_")
#' mod <- build_model_logistf(tmp, formula)
#' tmp <- mod$get_coefficients()
#' stopifnot(nrow(tmp) == 30)
#' mod$coef_histogram()
#' mod$coef_pairs()
#' mod$get_anova()
#' mod$coef_volcano()
#' mod$anova_histogram()
#' mod$write_coef_figures(tempdir())

ModelFirth <- R6::R6Class(
  "ModelFirth",
  inherit = ModelInterface,
  public = list(
    #' @field modelDF data.frame with modelling data and model.
    models = NULL,
    #' @field modelName name of model
    modelName = character(),
    #' @field subject_Id e.g. protein_Id
    subject_Id = character(),
    #' @field model_strategy function to create the models
    #' @field anova_df function to compute anova
    anova_df = NULL,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param modelDF dataframe with modelling results
    #' @param model_strategy model_strategy see \code{\link{strategy_lmer}}
    #' @param modelName name of model
    #' @param subject_Id subject column name
    #' @param p.adjust method to adjust p-values
    #'
    initialize = function(models,
                          modelName = "modelFirth",
                          subject_Id = "protein_Id",
                          p.adjust = prolfqua::adjust_p_values){
      self$models = models
      self$modelName = modelName
      self$subject_Id = subject_Id
      self$p.adjust = p.adjust

    },
    #' @description
    #' return model coefficient table
    get_coefficients = function(){
      lmermodel <- "linear_model"
      # Extract coefficients
      .coef_df <-  function(x){
        object <- x
        chi2 <- qchisq(1 - object$prob, 1)
        out <- cbind(object$coefficients, diag(object$var)^0.5,
                     object$ci.lower,
                     object$ci.upper,
                     chi2,
                     object$prob,
                     ifelse(object$method.ci=="Wald", 1, ifelse(object$method.ci=="-", 3, 2)))
        dimnames(out) <- list(names(object$coefficients), c("Estimate", "se(coef)", paste(c("lower", "upper"), 1 - object$alpha), "Chisq", "p", "method"))
        out <- data.frame(factor = row.names(out), out);
        return(out)
      }


      models <- dplyr::bind_rows(self$models[[1]]$modelDF, self$models[[2]]$modelDF)
      res <- vector(mode = "list", nrow(models))
      #for(i in 1:nrow(models)){
      #  print(i)
      #  x <- models$linear_model[[i]]
      #  res[[i]] <- .coef_df(x)
      #}
      Model_Coeff <- models |>
        dplyr::mutate(!!"Coeffs_model" := purrr::map( !!sym(lmermodel),  .coef_df ))
      Model_Coeff <- Model_Coeff |>
        dplyr::select(!!!syms(self$subject_Id), !!sym("Coeffs_model"), isSingular, nrcoef)
      Model_Coeff <- tidyr::unnest_legacy(Model_Coeff)
      if(!is.null(self$models$hkey)){
        Model_Coeff <- Model_Coeff |> dplyr::filter(!grepl(self$models$hkey, factor))
      }
      return(Model_Coeff)
    },
    #' @description
    #' return anova table
    get_anova = function(){
      warning("method not implemented!")
      return(NULL)
    },

    #' @description
    #' writes model coefficients to file
    #' @param path folder to write to
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write_coefficients  = function(path, format = "xlsx"){
      lfq_write_table(self$get_coefficients(),
                      path = path,
                      name  = paste0("Coef_",self$modelName),
                      format = format)
    },
    #' @description
    #' histogram of model coefficient
    coef_histogram = function(){
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      ## Coef_Histogram
      fname_histogram_coeff_p.values <- paste0("Coef_Histogram_",self$modelName,".pdf")
      histogram_coeff_p.values <- ggplot(data = Model_Coeff, aes(x = p, group = factor)) +
        geom_histogram(breaks = seq(0,1,by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = histogram_coeff_p.values, name = fname_histogram_coeff_p.values))
    },
    #' @description
    #' volcano plot of non intercept coefficients
    coef_volcano = function(){
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      fname_VolcanoPlot <- paste0("Coef_VolcanoPlot_",self$modelName,".pdf")
      VolcanoPlot <- Model_Coeff |>
        dplyr::filter(factor != "(Intercept)") |>
        prolfqua::multigroup_volcano(
          effect = "Estimate",
          significance = "p",
          contrast = "factor",
          label = "subject_Id" ,
          xintercept = c(-1, 1) ,
          colour = "isSingular" )
      return(list(plot = VolcanoPlot, name = fname_VolcanoPlot))
    },
    #' @description
    #' pairs-plot of coefficients
    coef_pairs = function(){
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      ## Coef_Pairsplot
      forPairs <- Model_Coeff |>
        dplyr::select(all_of(c("subject_Id" , "factor" ,  "Estimate") )) |>
        tidyr::pivot_wider(names_from = "factor", values_from = "Estimate" )
      fname_Pairsplot_Coef <- paste0("Coef_Pairsplot_", self$modelName,".pdf")
      #Pairsplot_Coef <-  GGally::ggpairs(forPairs, columns = 2:ncol(forPairs))
      return(list(plot = forPairs, name = fname_Pairsplot_Coef))

    },
    #' @description
    #' histogram of ANOVA results
    #' @param what show either "Pr..F." or "FDR.Pr..F."
    anova_histogram = function(what=c("p.value", "FDR")){
      warning("not implemented")
      return(NULL)
    },
    #' @description
    #' write figures related to ANOVA into pdf file
    #' @param path folder name
    #' @param width figure width
    #' @param height figure height
    #'
    write_anova_figures = function(path, width = 10, height =10){
      warning("not implemented")
      return(NULL)
    },
    #' @description
    #' write figures related to Coefficients into pdf file
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
