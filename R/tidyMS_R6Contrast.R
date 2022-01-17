# ContrastsInterface ----
#' ContrastsInterface
#' @export
ContrastsInterface <- R6::R6Class(
  "ContrastsInterface",
  public = list(
    #' @description
    #' get table with sides of the contrast
    get_contrast_sides = function(){stop("get_contrast_sides not implmented")},
    #' @description
    #' get table with contrast results (similar to limma topTable function)
    get_contrasts = function(){stop("get_contrasts not implmented")},
    #' @description
    #' initialize plotter
    get_Plotter = function(){stop("get_Plotter not implmented.") },
    #' @description
    #' create wide representation of data.
    to_wide = function(){stop("to_wide not implemented.")}
  )
)


.requiredContrastColumns <- c("contrast" , "c1" , "c2" ,
                              "c1_name" , "c2_name" , "sigma", "df",
                              "diff", "statistic", "p.value",
                              "conf.low", "conf.high", "FDR")

# ContrastsSimpleImpute----
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
#'
#' bb <- prolfqua_data('data_ionstar')$normalized()
#' configur <- bb$config$clone(deep=TRUE)
#' configur$table$hierarchyDepth <- 2
#' data <- bb$data
#' lfqdata <- LFQData$new(data, configur)
#' Contrasts <- c("dilution.b-a" = "dilution.b - dilution.a",
#' "dilution.c-e" = "dilution.c - dilution.b")
#' #ContrastsSimpleImpute$undebug("get_contrasts")
#' tmp <- ContrastsSimpleImpute$new(lfqdata, contrasts = Contrasts)
#' bb <- tmp$get_contrasts()
#' tmp$get_contrast_sides()
#' tmp$to_wide()
#' pl <- tmp$get_Plotter()
#' pl$histogram()
#' pl$histogram_diff()
#' pl$ma_plot()
#' pl$volcano()
#'
#' tmp <- ContrastsSimpleImpute$new(lfqdata, contrasts = Contrasts, method = "V2")
#' pl <- tmp$get_Plotter()
#' pl$histogram()
#'
ContrastsSimpleImpute <- R6::R6Class(
  "ContrastsSimpleImpute",
  inherit = ContrastsInterface,
  private = list(
    method = "V1"
  ),
  public = list(
    #' @field subject_Id subject_id e.g. protein_ID column
    subject_Id = character(),
    #' @field contrasts array with contrasts (see example)
    contrasts = character(),
    #' @field modelName model name
    modelName = character(),
    #' @field contrast_result data frame with results of contrast computation
    contrast_result = NULL,
    #' @field lfqdata data frame
    lfqdata = NULL,
    #' @field confint confidence interval
    confint = 0.95,
    #' @field p.adjust funciton to adjust p-values
    p.adjust = NULL,
    #' @field global Take global or local values for imputation
    global = logical(),
    #' @field probs qunatile to estimate missing values from.
    probs = 0.03,
    #' @description
    #' initialize
    #' @param lfqdata LFQData
    #' @param contrasts array of contrasts (see example)
    #' @param confint confidence interval
    #' @param p.adjust method for p-value adjustement - default benjamini hochberg
    #' @param modelName default "groupAverage"
    #' @param method internal default V1
    #' @param probs which quantile to use for imputation default 0.3
    #' @param global default TRUE use all or per condition data to impute from
    initialize = function(lfqdata,
                          contrasts,
                          confint = 0.95,
                          p.adjust = prolfqua::adjust_p_values,
                          modelName = "groupAverage",
                          method = "V1",
                          probs = 0.03,
                          global = TRUE){
      self$subject_Id = lfqdata$config$table$hkeysDepth()
      self$contrasts = contrasts
      self$modelName = modelName
      self$lfqdata = lfqdata
      self$confint = confint
      self$p.adjust = p.adjust
      private$method = method
      self$global  = global
      self$probs = probs
    },
    #' @description
    #' get contrasts sides
    #'
    get_contrast_sides = function(){
      if (is.null(self$contrast_result)) {
        self$get_contrasts()
      }
      self$contrast_result |> dplyr::select(contrast,c1 = c1_name,c2 = c2_name) |> distinct()
    },
    #' @description
    #' table with results of contrast computation
    #' @param all FALSE, do not show all columns (default)
    get_contrasts = function(all = FALSE){
      if (is.null(self$contrast_result)) {
        result = get_imputed_contrasts(self$lfqdata$data, self$lfqdata$config, self$contrasts, probs = self$probs, global = self$global)

        if (self$lfqdata$config$table$hierarchyDepth < length(self$lfqdata$config$table$hierarchyKeys())) {
          stop("hierarchy depth < hierarchyKeys(). Please aggregate first.")
        } else {
          # compute statistics using pooled variance
          result$isSingular <- TRUE
          #result = get_imputed_contrasts(transformed$data, transformed$config, Contrasts)
          result <- select(result , -all_of(c("n","estimate_mad")))

          var = summarize_stats(self$lfqdata$data, self$lfqdata$config)
          pooled <- poolvar(var, self$lfqdata$config, method = self$method)
          pooled <- dplyr::select(pooled ,-all_of(c(self$lfqdata$config$table$fkeysDepth()[1],"var")))

          result <- dplyr::inner_join(result, pooled, by = self$lfqdata$config$table$hkeysDepth())
          result <- dplyr::mutate(result, statistic = .data$estimate_median / .data$sdT,
                                  p.value = 2*pt(abs(.data$statistic), df = .data$df, lower.tail = FALSE))

          prqt <- -qt((1 - self$confint)/2, df = result$df)
          result$conf.low <- result$estimate_median  - prqt * (result$sdT)
          result$conf.high <- result$estimate_median + prqt * (result$sdT)
          result <- self$p.adjust(result, column = "p.value", group_by_col = "contrast", newname = "FDR")
          if (!all) {
            result <- select(result, -all_of( c("isSingular", "not_na" , "mean" ,"n.groups", "n", "meanAll", "sdT") ) )
          }

        }

        result <- result |> rename(diff = estimate_median, sigma = sd)
        result <- mutate(result, modelName = self$modelName, .before = 1)
        self$contrast_result <- result
      }
      res <- ungroup(self$contrast_result)
      stopifnot(all(.requiredContrastColumns %in% colnames(res)))
      invisible(res)
    },
    #' @description
    #' get Contrasts_Plotter
    #' @return Contrast_Plotter
    get_Plotter = function(){
      res <- Contrasts_Plotter$new(
        self$get_contrasts(),
        subject_Id = self$subject_Id,
        volcano = list(list(score = "FDR", thresh = 0.1)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "diff",
        contrast = "contrast")
      return(res)
    },
    #' @description
    #' convert contrast results to wide format
    #' @param columns value column default p.value
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("diff", columns),
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
#'
#'
#' istar <- prolfqua_data('data_ionstar')$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#' strategy_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id) + (1 | sampleName)")
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
#'     "dil.e_vs_b" = "dilution.e - dilution.b" )
#'  #prolfqua::Contrasts$debug("get_linfct")
#' #prolfqua::Contrasts$debug("get_contrasts")
#' contrastX <- prolfqua::Contrasts$new(mod,
#'  Contr)
#'
#' contrastX$get_contrast_sides()
#'
#' contrastX$get_linfct()
#' xx <- contrastX$get_contrasts()
#'
#'
#' contrastX$get_contrasts()
#' contrastX$to_wide()
#'
Contrasts <- R6::R6Class(
  "Contrast",
  inherit = ContrastsInterface,
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
    #' @field p.adjust function to adjust p-values (default prolfqua::adjust_p_values)
    p.adjust = NULL,
    #' @field contrast_result data frame containing results of contrast computation
    contrast_result = NULL,
    #' @field global use a global linear function (determined by get_linfct)
    global = TRUE,
    #' @description
    #' initialize
    #' create Contrast
    #' @param model a dataframe with a structure similar to that generated by \code{\link{build_model}}
    #' @param contrasts a character vector with contrast specificiation
    #' @param p.adjust function to adjust the p-values
    #' @param global development/internal argument (if FALSE determine linfct for each model.)
    #' @param modelName name of contrast method, default WaldTest
    initialize = function(model,
                          contrasts,
                          p.adjust = prolfqua::adjust_p_values,
                          global = TRUE,
                          modelName = "WaldTest"
    ){
      self$models = model$modelDF
      self$contrasts = contrasts
      self$contrastfun = model$modelFunction$contrast_fun
      self$modelName =  modelName
      self$subject_Id = model$subject_Id
      self$p.adjust = p.adjust
      self$global = global
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function(){
      # extract contrast sides
      tt <- self$contrasts[grep("-",self$contrasts)]
      tt <- tibble(contrast = names(tt) , rhs = tt)
      tt <- tt |> mutate(rhs = gsub("[` ]","",rhs)) |>
        tidyr::separate(rhs, c("c1", "c2"), sep = "-")
      return(tt)
    },
    #' @description
    #' get linear functions from contrasts
    #' @param global logical TRUE - get the a linear functions for all models, FALSE - linear function for each model
    get_linfct = function(global = TRUE){
      linfct <- function(model, contrast){
        linfct <- linfct_from_model(model, as_list = FALSE)
        linfct <- unique(linfct) # needed for single factor models
        linfct_A <- linfct_matrix_contrasts(linfct, self$contrasts)
        return( rbind( linfct, linfct_A ) )
      }
      if (global) {
        models <- self$models |> dplyr::filter(exists_lmer == TRUE)
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
    #' @return data.frame with contrasts
    #'
    get_contrasts = function(all = FALSE){
      if (is.null(self$contrast_result) ) {

        message("determine linear functions:")
        linfct <- self$get_linfct(global = self$global)
        contrast_sides <- self$get_contrast_sides()
        message("compute contrasts:")
        contrast_result <- contrasts_linfct(self$models,
                                            linfct,
                                            subject_Id = self$subject_Id,
                                            contrastfun = self$contrastfun )

        contrast_result <- dplyr::rename(contrast_result, contrast = lhs, diff = estimate)

        xx <- dplyr::select(contrast_result, self$subject_Id, "contrast", "diff")
        xx <- xx |> tidyr::pivot_wider(names_from = "contrast", values_from = "diff")

        contrast_result <- dplyr::filter(contrast_result, contrast %in% names(self$contrasts))

        get_contrast_cols <- function(i, contrast_results , contrast_table , subject_ID ){
          data.frame(lhs = contrast_table[i, "contrast"],
                     dplyr::select_at(contrast_results, c( subject_ID, unlist(contrast_table[i,c("c1","c2")]))),
                     c1_name = contrast_table[i,"c1", drop = TRUE],
                     c2_name = contrast_table[i,"c2", drop = TRUE], stringsAsFactors = FALSE)
        }

        contrast_sides <- purrr::map_df(seq_len(nrow(contrast_sides)),
                                        get_contrast_cols,
                                        xx,
                                        contrast_sides,
                                        self$subject_Id)
        contrast_result <- inner_join(contrast_sides, contrast_result)

        contrast_result <- self$p.adjust(contrast_result,
                                         column = "p.value",
                                         group_by_col = "contrast",
                                         newname = "FDR")
        contrast_result <- mutate(contrast_result, modelName = self$modelName, .before = 1)
        self$contrast_result <- dplyr::ungroup(contrast_result )
      }

      res <- if (!all) {
        self$contrast_result |>
          select( -all_of(c("sigma.model",
                            "df.residual.model",
                            "std.error",
                            "isSingular")))

      }else{
        self$contrast_result
      }
      res <- ungroup(res)

      stopifnot(all(.requiredContrastColumns %in% colnames(res)))
      return(res)

    },
    #' @description
    #' return \code{\link{Contrasts_Plotter}}
    #' creates Contrast_Plotter
    #' @param FCthreshold fold change threshold to show in plots
    #' @param FDRthreshold FDR threshold to show in plots
    #' @return \code{\link{Contrasts_Plotter}}
    get_Plotter = function(
      FCthreshold = 1,
      FDRthreshold = 0.1
    ){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(list(score = "FDR", thresh = FDRthreshold)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        score = list(list(score = "statistic", thresh = 5)),
        modelName = self$modelName,
        estimate = "diff",
        contrast = "contrast"
        )
      return(res)
    },
    #' @description
    #' convert to wide format
    #' @param columns value column default p.value
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("diff", columns),
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
#' istar <- prolfqua_data('data_ionstar')$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   strategy_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id) + (1|sampleName)")
#' pepIntensity <- istar_data
#' config <- istar$config$clone(deep = TRUE)
#' config$table$hkeysDepth()
#' mod <- build_model(
#'  pepIntensity,
#'  modelFunction,
#'  subject_Id = config$table$hkeysDepth())
#'
#'  Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")
#'  contrast <- prolfqua::Contrasts$new(mod,
#'  Contr)
#'  #ContrastsModerated$debug("get_Plotter")
#'  contrast <- ContrastsModerated$new(contrast)
#'  bb <- contrast$get_contrasts()
#'
#'  plotter <- contrast$get_Plotter()
#'  plotter$histogram()
#'  plotter$ma_plot()
#'  plotter$volcano()
#'
#' #bb |> dplyr::rename(log2FC = estimate, mean_c1 = c1, mean_c2 = c2)
ContrastsModerated <- R6::R6Class(
  classname = "ContrastsModerated",
  inherit = ContrastsInterface,
  public = list(
    #' @field Contrast Class implementing the Contrast interface
    Contrast = NULL,
    #' @field modelName name of model
    modelName = character(),
    #' @field subject_Id columns with subject_Id (proteinID)
    subject_Id = character(),
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param Contrast class implementing the ContrastInterface
    #' @param modelName name of the model
    #' @param p.adjust function to adjust p-values - default BH
    initialize = function(Contrast,
                          modelName = paste0(Contrast$modelName, "_moderated"),
                          p.adjust = prolfqua::adjust_p_values
    ){
      self$Contrast = Contrast
      self$subject_Id = Contrast$subject_Id
      self$modelName = modelName
      self$p.adjust = p.adjust
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function(){
      self$Contrast$get_contrast_sides()
    },
    #' @description
    #' get linear functions from contrasts
    #' @param global logical TRUE - get the a linear functions for all models, FALSE - linear function for each model
    get_linfct = function(global = TRUE){
      self$Contrast$get_linfct()
    },
    #' @description
    #' applies limma moderation
    #' @seealso \code{\link{moderated_p_limma_long}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      contrast_result <- self$Contrast$get_contrasts(all = FALSE)
      contrast_result <- moderated_p_limma_long(
        contrast_result ,
        group_by_col = "contrast",
        estimate = "diff")
      if (!all) {
        contrast_result <- contrast_result |> select(-c( "sigma","df",
                                                          "statistic", "p.value","conf.low","conf.high",
                                                          "FDR",  "moderated.df.prior" ,
                                                          "moderated.var.prior"))
        contrast_result <- contrast_result |> mutate(sigma = sqrt(moderated.var.post),.keep = "unused")
        contrast_result <- contrast_result |> rename(
          conf.low = "moderated.conf.low",
          conf.high = "moderated.conf.high",
          statistic = "moderated.statistic" ,
          df = "moderated.df.total",
          p.value = "moderated.p.value"
        )
        contrast_result <- self$p.adjust(contrast_result, column = "p.value",
                                         group_by_col = "contrast",
                                         newname = "FDR")
      }else{
        contrast_result <- self$p.adjust(contrast_result,
                                         column = "moderated.p.value",
                                         group_by_col = "contrast",
                                         newname = "FDR.moderated")
      }
      contrast_result <- mutate(contrast_result,modelName = self$modelName, .before  = 1)
      contrast_result <- dplyr::ungroup(contrast_result)
      stopifnot(all(.requiredContrastColumns %in% colnames(contrast_result)))

      return(contrast_result)
    },
    #' @description
    #' get \code{\link{Contrasts_Plotter}}
    #' @param FCthreshold fold change threshold to show in plots
    #' @param FDRthreshold FDR threshold to show in plots
    #'
    get_Plotter = function(
      FCthreshold = 1,
      FDRthreshold = 0.1
    ){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(list(score = "FDR", thresh = FDRthreshold)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        score = list(list(score = "statistic", thresh = 5)),
        modelName = self$modelName,
        estimate = "diff",
        contrast = "contrast"
      )
      return(res)
    },
    #' @description
    #' convert to wide format
    #' @param columns value column default moderated.p.value
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("diff", columns),
                                                     contrast = 'contrast')
      return(contrasts_wide)
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
#' istar <- prolfqua_data('data_ionstar')$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelFunction <-
#'   strategy_lm("transformedIntensity  ~ dilution.")
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
#'  Contr <- c("dil.b_vs_a" = "dilution.b - dilution.a")
#'
#'
#'  contr <- prolfqua::Contrasts$new(mod, Contr)
#'  dim(contr$get_contrasts())
#'  contrM <- prolfqua::ContrastsModerated$new(contr)
#'  dim(contrM$get_contrasts())
#'  contrast <- prolfqua::ContrastsROPECA$new(contrM)
#'  contrast$get_contrasts()
#'  #ContrastsROPECA$debug("to_wide")
#'  contrast <- prolfqua::ContrastsROPECA$new(contr)
#'  tmp <- contrast$get_contrasts()
#'  dim(tmp)
#'  pl <- contrast$get_Plotter()
#'  contrast$to_wide()
#'  pl$histogram()
#'  pl$ma_plot()
#'
ContrastsROPECA <- R6::R6Class(
  "ContrastsROPECA",
  inherit = ContrastsInterface,
  public = list(
    #' @field Contrast Contrast
    Contrast = NULL,
    #' @field contrast_result contrast result
    contrast_result = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id columns with protein ID's
    subject_Id = character(),
    #' @field p.adjust method to use for p.value adjustment
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param Contrast e.g. instance of Contrasts class, or ContrastsModerated
    #' @param modelName default ROPECA
    #' @param p.adjust function to use for p.value adjustement
    initialize = function(Contrast,
                          modelName = "ROPECA",
                          p.adjust = prolfqua::adjust_p_values
    ){
      self$Contrast = Contrast
      stopifnot(length(Contrast$subject_Id) > 1)
      self$modelName = modelName
      self$subject_Id = Contrast$subject_Id
      self$p.adjust = p.adjust
    },
    #' @description
    #' show names of contrasts
    #' @return data.frame
    get_contrast_sides = function(){
      if (is.null(self$contrast_result)) {
        self$get_contrasts()
      }
      self$contrast_result |> dplyr::select(contrast,c1 = c1_name,c2 = c2_name) |> distinct()
    },
    #' @description
    #' get linear function used to determine contrasts
    #' @return data.frame
    get_linfct = function(){
      self$contrast$get_linfct()
    },
    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    #' @return data.frame
    get_contrasts = function(all = FALSE){
      if (is.null(self$contrast_result)) {
        contrast_result <- self$Contrast$get_contrasts(all = FALSE)
        contrast_result <- summary_ROPECA_median_p.scaled(
          contrast_result,
          contrast = "contrast",
          subject_Id = self$subject_Id[length(self$subject_Id) - 1],
          estimate = "diff",
          statistic = "statistic",
          p.value = "p.value",
          max.n = 10)
        contrast_result <- dplyr::rename(contrast_result, diff = "estimate")
        contrast_result <- self$p.adjust(contrast_result,
                                         column = "beta.based.significance",
                                         group_by_col = "contrast",
                                         newname = "FDR.beta.based.significance")
        contrast_result <- self$p.adjust(contrast_result,
                                         column = "median.p.value",
                                         group_by_col = "contrast",
                                         newname = "FDR.median.p.value")
        contrast_result <- mutate(contrast_result,modelName = self$modelName, .before  = 1)
        self$contrast_result <- contrast_result

      }

      if (!all) {
        res <- select(self$contrast_result ,
                      -all_of(c( "n_not_na", "mad.estimate",
                                 "n.beta", "isSingular",
                                 "median.p.scaled","median.p.value",
                                 "FDR.median.p.value")) )
      }else{
        res <- self$contrast_result
      }

      return(res)
    },
    #' @description
    #' get \code{\link{Contrasts_Plotter}}
    #' @return \code{\link{Contrasts_Plotter}}
    #' @param FDRthreshold FDR threshold
    #' @param FCthreshold FC threshold
    get_Plotter = function(FDRthreshold = 0.1,
                           FCthreshold = 2){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id[1],
        fcthresh = FCthreshold,
        volcano = list(list(score = "FDR.beta.based.significance", thresh = FDRthreshold)),
        histogram = list(list(score = "beta.based.significance", xlim = c(0,1,0.05)),
                         list(score = "FDR.beta.based.significance", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "diff",
        contrast = "contrast")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default beta.based.significance
    #' @return data.frame
    to_wide = function(columns = c("beta.based.significance", "FDR.beta.based.significance")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(
        contrast_minimal,
        subject_Id = self$subject_Id[length(self$subject_Id) - 1],
        columns = c("diff", columns),
        contrast = 'contrast')
      return(contrasts_wide)
    }

  ))


# ContrastsSaintExpress ----
#' ContrastsSaintExpress Wrapper to results produced by SaintExpress (list.txt file)
#'
#' @export
#' @family modelling
#'
ContrastsSaintExpress <- R6::R6Class(
  "ContrastsSaint",
  inherit = ContrastsInterface,
  public = list(
    #' @field contrast_result data.frame with the contrast computation results
    contrast_result = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id subject id defualt 'Prey'
    subject_Id = character(),
    #' @description
    #' initialize
    #' @param contrastsdf return value of \code{\link{runSaint}}
    #' @param subject_Id default "Prey"
    #' @param modelName name of model
    initialize = function(contrastsdf,
                          subject_Id = "Prey",
                          modelName = "ContrastSaint"
    ){


      self$contrast_result = contrastsdf
      self$subject_Id = subject_Id
      self$modelName = modelName

      if ( "AvgIntensity" %in% colnames(contrastsdf)) {
        self$contrast_result <- contrastsdf |> mutate(log2FC = log2(FoldChange),
                                                       c1_name = "Control",
                                                       c2_name = Bait,
                                                       c1 = log2(AvgIntensity) - log2(FoldChange)/2,
                                                       c2 = log2(AvgIntensity) + log2(FoldChange)/2 ,
                                                       modelName = modelName)

      }else{
        self$contrast_result <- contrastsdf |> mutate(log2FC = log2(FoldChange),
                                                       c1_name = "Control",
                                                       c2_name = Bait,
                                                       c1 = log2(AvgSpec) - log2(FoldChange)/2,
                                                       c2 = log2(AvgSpec) + log2(FoldChange)/2 ,
                                                       modelName = modelName)

      }

    },
    #' @description
    #' show contrasts
    #' @return data.frame
    get_contrast_sides = function(){
      self$contrast_result |> dplyr::select(Bait,c1 = c1_name,c2 = c2_name) |> distinct()
    },
    #' @description
    #' no available for SaintExpress
    #'
    get_linfct = function(){
      NULL
    },

    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      res <- self$contrast_result |> select(
        all_of(c(self$subject_Id,
                 "modelName",
                 "Bait",
                 "c1_name",
                 "c2_name",
                 "c1",
                 "c2",
                 "log2FC",
                 "SaintScore",
                 "BFDR"
        )))
      res
    },
    #' @description get \code{\link{Contrasts_Plotter}}
    #' @param FCthreshold fold change threshold to show
    #' @param SaintScore SaintScore threshold to show in the heatmap.
    #' @param BFDRthreshold BDRF threshold
    #' @return \code{\link{Contrasts_Plotter}}
    get_Plotter = function(FCthreshold = 1, SaintScore = 0.75, BFDRthreshold = 0.1){
      res <- Contrasts_Plotter$new(
        self$contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(list(score = "BFDR", thresh = BFDRthreshold)),
        histogram = list(list(score = "BFDR", xlim = c(0,1,0.05)), list(score = "SaintScore", xlim = c(0,1,0.05))),
        score = list(list(score = "SaintScore", thresh = SaintScore )),
        modelName = self$modelName,
        estimate = "log2FC",
        contrast = "Bait")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default SaintScore, BFDR
    #' @return data.frame
    to_wide = function(columns = c("SaintScore", "BFDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("log2FC", columns),
                                                     contrast = 'Bait')
      return(contrasts_wide)
    }
  ))



# ContrastsTable -----

#'
#' ContrastTable (place holder future baseclass?)
#'
#'
#' @export
#' @family modelling
#' @examples
#'
#' bb <- prolfqua_data('data_ionstar')$normalized()
#' configur <- bb$config$clone(deep=TRUE)
#' configur$table$hierarchyDepth <- 2
#' data <- bb$data
#' lfqdata <- LFQData$new(data, configur)
#' Contrasts <- c("dilution.b-a" = "dilution.b - dilution.a",
#' "dilution.c-e" = "dilution.c - dilution.b")
#' tmp <- ContrastsSimpleImpute$new(lfqdata, contrasts = Contrasts)
#' ctr <- tmp$get_contrasts()
#' xcx <- ContrastsTable$new(ctr, subject_Id = tmp$subject_Id, modelName = tmp$modelName)
#' xcx$get_Plotter()$volcano()
#'
ContrastsTable <- R6::R6Class(
  "ContrastsTable",
  inherit = ContrastsInterface,
  public = list(
    #' @field contrast_result contrast results
    contrast_result = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id default protein_Id
    subject_Id = character(),
    #' @description
    #' intitialize
    #' @param contrastsdf data.frame
    #' @param subject_Id default protein_Id
    #' @param modelName default ContrastTable
    initialize = function(contrastsdf,
                          subject_Id = "protein_Id",
                          modelName = "ContrastTable"
    ){
      self$contrast_result = contrastsdf
      self$subject_Id = subject_Id
      self$modelName = modelName
    },
    #' @description
    #' return sides of contrast
    #' @return data.frame
    get_contrast_sides = function(){
      self$contrast_result |> dplyr::select(contrast,c1 = c1_name,c2 = c2_name) |> distinct()
    },
    #' @description not implemented
    get_linfct = function(){
      NULL
    },

    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      self$contrast_result
    },
    #' @description
    #' get \code{\link{Contrasts_Plotter}}
    #' @param FCthreshold fold change threshold
    #' @param FDRthreshold fdr threshold
    #' @return \code{\link{Contrasts_Plotter}}
    #'
    get_Plotter = function(FCthreshold = 1, FDRthreshold = 0.1){
      res <- Contrasts_Plotter$new(
        self$contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(list(score = "FDR", thresh = FDRthreshold)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "diff",
        contrast = "contrast")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default beta.based.significance
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("diff", columns),
                                                     contrast = 'contrast')
      return(contrasts_wide)
    }
  )
)




# Contrasts_Plotter ----
#' plot contrasts
#' @export
#' @family modelling
#' @family plotting
#' @examples
#'
#'
#'
#' istar <- prolfqua_data('data_ionstar')$normalized()
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 100))
#' modelName <- "Model"
#' modelFunction <-
#'   strategy_lmer("transformedIntensity  ~ dilution. + (1 | peptide_Id)",
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
#'   "dil.e_vs_a" = "dilution.e - dilution.a"
#'   ,"dil.e_vs_b" = "dilution.e - dilution.b",
#'   "dil.c_vs_b" = "dilution.c - dilution.b"
#'  )
#' #Contrasts$debug("get_contrasts")
#' contrast <- prolfqua::Contrasts$new(mod,
#'   Contr)
#' tmp <- contrast$get_contrasts()
#' #Contrasts_Plotter$debug("score_plot")
#'
#'
#' cp <- Contrasts_Plotter$new(tmp ,
#'  contrast$subject_Id,
#' volcano = list(list(score = "FDR", thresh = 0.1)),
#' histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
#'                  list(score = "FDR", xlim = c(0,1,0.05))),
#' score =list(list(score = "statistic",  fcthresh = 2, thresh = 5)))
#' p <- cp$score_plot(legend=FALSE)
#' cp$score_plotly()
#' p <- cp$histogram()
#' p <- cp$histogram_diff()
#' res <- cp$volcano()
#' length(res)
#' res
#' respltly <- cp$volcano_plotly()
#' length(respltly)
#' cp$ma_plot()
#' cp$ma_plotly()
#' res  <- cp$barplot_threshold()
Contrasts_Plotter <- R6::R6Class(
  "Contrasts_Plotter",
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
    #' @field volcano_spec volcano plot specification
    volcano_spec = NULL,
    #' @field score_spec score plot specification
    score_spec = NULL,
    #' @field histogram_spec plot specification
    histogram_spec = NULL,
    #' @field fcthresh fold change threshold
    fcthresh = 1,
    #' @description
    #' create Crontrast_Plotter
    #' @param contrastDF frame with contrast data
    #' @param subject_Id columns containing subject Identifier
    #' @param volcano which score to plot and which ablines to add.
    #' @param histogram which scores to plot and which range (x) should be shown.
    #' @param score score parameters
    #' @param fcthresh default 1 (log2 FC threshold)
    #' @param modelName name of model
    #' @param estimate fold change (difference) estimate column
    #' @param contrast contrast column
    initialize = function(contrastDF,
                          subject_Id,
                          volcano = list(list(score = "FDR", thresh = 0.1)),
                          histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                                           list(score = "FDR", xlim = c(0,1,0.05))),
                          score = list(list(score = "statistic",  thresh = NULL)),
                          fcthresh = 1,
                          modelName = "Model",
                          estimate = "diff",
                          contrast = "contrast"
    ){
      self$contrastDF <- tidyr::unite(
        contrastDF,
        "subject_Id", subject_Id, sep = "~", remove = FALSE)

      self$modelName  = modelName
      self$subject_Id = subject_Id
      self$estimate = estimate
      self$volcano_spec = volcano
      self$score_spec = score
      self$histogram_spec = histogram
      self$fcthresh = fcthresh
      self$contrast = contrast
    },
    #' @description
    #' plot histogram of selected scores (e.g. p-value, FDR, t-statistics)
    histogram = function(){
      res <- list()
      if (length(self$histogram_spec) > 0) {
        for (spec in self$histogram_spec) {
          fig <- private$.histogram(score = spec )
          res[[spec$score]] <- fig
        }
        return(res)
      }
    },
    #' @description
    #' plot histogram of estimated fold change
    #' @param binwidth with of bin in histogram
    histogram_estimate = function(binwidth = 0.05){
      re <- range(self$contrastDF[[self$estimate]], na.rm = TRUE)
      re[1] <- floor(re[1])
      re[2] <- ceiling(re[2])
      fig <- private$.histogram(score = list(score =  self$estimate, xlim = c(re,binwidth)))
      return(fig)
    },
    #' @description
    #' plot histogram of differences (estimates) fold change
    #' @param binwidth with of bin in histogram
    histogram_diff = function(binwidth = 0.05){
      self$histogram_estimate(binwidth = binwidth)
    },
    #' @description
    #' volcano plots (fold change vs FDR)
    #' @param colour column name with color information default modelName
    volcano = function(colour = "modelName"){
      fig <- private$.volcano(self$contrastDF, self$volcano_spec, colour = colour )
      return(fig)
    },
    #' @description
    #' plotly volcano plots
    #' @param colour column in contrast matrix with colour coding
    #' @return list of ggplots
    volcano_plotly = function(colour = "modelName"){
      contrastDF <- self$contrastDF |> plotly::highlight_key(~ subject_Id)
      res <- private$.volcano(contrastDF, self$volcano_spec, colour = colour )
      for (i in seq_along(res)) {
        res[[i]] <- res[[i]] |> plotly::ggplotly(tooltip = "subject_Id")
      }
      return(res)
    },
    #' @description
    #' ma plot
    #' @param fc fold change abline
    #' @param colour column in contrast matrix with colour coding
    #' @param legend enable legend default TRUE
    #' @return ggplot
    ma_plot = function(fc, colour = "modelName", legend = TRUE){
      if ( missing(fc))
        fc <- self$fcthresh
      contrastDF <- self$contrastDF
      if (!is.null(contrastDF$c1) && !is.null(contrastDF$c2)) {
        # pdf version
        fig <- private$.ma_plot(contrastDF ,self$contrast, fc, colour = colour, legend = legend)
      }else{
        warning("no c1 c2 columns can't generate MA")
        fig <- NULL
      }
      return(fig)
    },

    #' @description
    #' ma plotly
    #' @param fc horizontal lines
    #' @param colour column in contrast matrix with colour coding.
    #' @param legend enable legend default TRUE
    #' @return list of ggplots
    ma_plotly = function(fc, colour = "modelName", legend = TRUE){
      # html version
      if (missing(fc))
        fc <- self$fcthresh
      contrastDF  <- self$contrastDF
      if (!is.null(contrastDF$c1) && !is.null(contrastDF$c2)) {
        contrastDF  <- contrastDF |>
          plotly::highlight_key(~subject_Id)
        fig_plotly <- private$.ma_plot(contrastDF, self$contrast, fc, colour = colour, legend = legend) |>
          plotly::ggplotly(tooltip = "subject_Id")

        return(fig_plotly)
      }else{
        return(NULL)
      }
    },
    #' @description
    #' plot a score against the log2 fc e.g. t-statistic
    #' @param scorespec list(score="statistics", fcthres = 2, thresh = 5)
    #' @param colour column with colour coding
    #' @param legend enable legend default TRUE
    #' @return list of ggplots
    score_plot = function(scorespec,  colour = "modelName", legend = TRUE ){
      if (!missing(scorespec)) {
        self$score_spec[[scorespec$score]] <- scorespec
      }
      res <- list()
      if (length(self$score_spec) > 0) {
        res <- private$.score_plot(
          self$contrastDF,
          self$score_spec,
          colour = colour,
          legend = legend )
      }
      return(res)
    },
    #' @description
    #' plot a score against the log2 fc e.g. t-statistic
    #' @param scorespec list(score="statistics", fcthres = 2, thresh = 5)
    #' @param colour column with colour coding
    #' @param legend enable legend default TRUE
    #' @return list of ggplots
    score_plotly = function(scorespec,  colour = "modelName", legend = TRUE ){
      if (!missing(scorespec)) {
        self$score_spec[[scorespec$score]] <- scorespec
      }
      contrastDF <- self$contrastDF |> plotly::highlight_key(~ subject_Id)
      res <- private$.score_plot(
        contrastDF,
        self$score_spec,
        colour = colour,
        legend = legend )

      for (i in seq_along(res)) {
        res[[i]] <- plotly::ggplotly(res[[i]], tooltip = "subject_Id")
      }
      return(res)
    },
    #' @description
    #' shows the number of significant protens per contrasts
    #' @return list which contains ggplots and summary tables
    barplot_threshold = function(){
      resBar <- list()
      for (i in seq_along(self$volcano_spec) ) {
        scN <- self$volcano_spec[[i]]$score
        scT <- self$volcano_spec[[i]]$thresh
        filt <- self$contrastDF  |> filter(!is.na(!!sym(scN)) & !!sym(scN)  < scT)
        sumC <- filt |> group_by(contrast, modelName) |> dplyr::summarize(n = n())
        p <- ggplot(sumC, aes(x = contrast, y = n, fill = modelName)) + geom_bar(position = "stack", stat = "identity")
        resBar[[scN]] <- list(plot = p, summary = sumC)
      }
      return(resBar)
    }
  ),
  private = list(
    .volcano = function(contrasts, scores,  colour = NULL, legend = TRUE){
      fig <- list()
      for (score in scores) {
        column <- score$score
        p <- prolfqua:::.multigroupVolcano(
          contrasts,
          effect = self$estimate,
          p.value = column,
          condition = self$contrast,
          text = "subject_Id",
          xintercept = c(-self$fcthresh, self$fcthresh),
          yintercept = score$thresh,
          colour = colour,
          scales = "free_y")
        if (!legend) {
          p <- p + guides(colour = "none")
        }
        fig[[column]] <- p
      }
      return(fig)
    },
    .histogram  = function(score){
      xlim = score$xlim
      score = score$score
      plot <- self$contrastDF |> ggplot(aes(x = !!sym(score))) +
        geom_histogram(breaks = seq(from = xlim[1], to = xlim[2], by = xlim[3])) +
        facet_wrap(vars(!!sym(self$contrast))) +
        theme_light()
      return(plot)
    },
    .ma_plot = function(x, contrast, fc, colour = NULL, legend = TRUE){
      p <- ggplot(x , aes(x = (c1 + c2)/2,
                          y = !!sym(self$estimate),
                          text = !!sym("subject_Id"),
                          colour = !!sym(colour))) +
        geom_point(alpha = 0.5) +
        scale_colour_manual(values = c("black", "green")) +
        facet_wrap(vars(!!sym(contrast))) + theme_light() +
        geom_hline(yintercept = c(-fc, fc), linetype = "dashed",colour = "red") +
        ylab("log fold change (M)") + xlab("mean log intensities (A)") +
        theme_light()
      if (!legend) {
        p <- p + guides(colour = "none")
      }
      return(p)
    },
    .score_plot = function(x, scores, colour, legend = TRUE){
      fig <- list()
      for (score in scores) {
        xlim = self$fcthresh
        ylim = score$thresh
        score = score$score
        scoreVal <- if ("data.frame" %in% class(x)) {
          x[[score]]
        } else {
          x$data()[[score]]
        }

        ylims <- c( sign(min(scoreVal, na.rm = TRUE)) * ylim, sign( max(scoreVal, na.rm = TRUE)) * ylim)
        p <- ggplot(x, aes(x = !!sym(self$estimate),
                           y = !!sym(score),
                           text = !!sym("subject_Id"),
                           colour = !!sym(colour))) +
          scale_colour_manual(values = c("black", "green")) +
          geom_point(alpha = 0.5) +
          facet_wrap(vars(!!sym(self$contrast))) +
          geom_hline(yintercept = c(0), colour = 1) +
          geom_vline(xintercept = c(0), colour = 1 ) +
          geom_hline(yintercept = ylims, colour = 2, linetype = "dashed") +
          geom_vline(xintercept = c(-xlim, xlim), colour = 2, linetype = "dashed" ) +
          theme_light()
        if (!legend) {
          p <- p + guides(colour = "none")
        }
        fig[[score]] <- p
      }
      return(fig)
    }
  )
)



# Merge contrasts ----
#' add contrast results from two different functions. Tipically used with Contrast and Cotnrast simple imputed.
#'
#' @param prefer contrasts to use preferentially
#' @param add contrasts to add from if missing in prefer
#' @param modelName name of the merged model default "mergedModel"
#' @export
#' @family modelling
#'
addContrastResults <- function(prefer, add, modelName = "mergedModel"){
  cA <- prefer$get_contrasts()
  cB <- add$get_contrasts()

  stopifnot(length(setdiff(colnames(cA), colnames(cB))) == 0)

  cA <- dplyr::filter(cA,!is.na(.data$statistic))
  moreID <- setdiff(distinct(select(cB, c(prefer$subject_Id, "contrast"))),
                    distinct(select(cA, c(add$subject_Id, "contrast"))))
  more <- inner_join(moreID, cB )

  sameID <- select(cA, c(add$subject_Id, "contrast"))
  same <- inner_join(sameID, cB)

  merged <- bind_rows(cA, more)

  if (prefer$modelName == add$modelName) {
    prefermodelName <- paste0(prefer$modelName, "_prefer")
    addmodelName <- paste0(add$modelName, "_add")
    cA$modelName <- prefermodelName
    more$modelName <- addmodelName
  } else {
    prefermodelName <- prefer$modelName
    addmodelName <- add$modelName
  }

  merged$modelName <- factor(merged$modelName,
                             levels = c(levels(factor(cA$modelName)), addmodelName))

  merged <- ContrastsTable$new(merged,
                               subject_Id = prefer$subject_Id,
                               modelName = paste0(prefermodelName,"_",addmodelName))
  more <- ContrastsTable$new(more , subject_Id = prefer$subject_Id, modelName = addmodelName)
  same <-  ContrastsTable$new(same , subject_Id = prefer$subject_Id, modelName = addmodelName)
  return(list(merged = merged, more = more, same = same))
}


