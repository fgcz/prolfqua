# ContrastsInterface ----
#' ContrastsInterface
#' @export
ContrastsInterface <- R6::R6Class(
  "ContrastsInterface",
  public = list(
    get_contrast_sides = function(){stop("get_contrast_sides not implmented")},
    get_contrasts = function(){stop("get_contrasts not implmented")},
    get_Plotter = function(){stop("get_Plotter not implmented.") },
    to_wide = function(){stop("to_wide not implemented.")}
  )
)
# summarise_missing_contrasts
#' @examples
#' ttd <- ionstar_bench_preprocess(prolfqua::data_benchmarkExample)
#' x <- .summarise_missing_contrasts(ttd$data)
#' x2 <- as_tibble(x$summary)
#'
#debug(.summarise_missing_contrasts)
.summarise_missing_contrasts <- function(data,
                                         hierarchy = c("protein_Id"),
                                         contrast = "contrast",
                                         what = "statistic") {
  data <- tidyr::complete(
    data,
    tidyr::nesting(!!!syms(contrast)),
    tidyr::nesting(!!!syms(hierarchy))
  )

  xxA <- data |>
    group_by_at(hierarchy) |>
    summarize(n = n(), nr_na = sum(is.na(!!sym(what))))
  summary <- xxA |> group_by(.data$nr_na) |> summarize(n = n())

  colnames(summary) <- c("nr_missing", paste(hierarchy, collapse = "_"))
  return(list(summary = summary, nr_na = xxA))
}

.requiredContrastColumns <- c("contrast" , "c1" , "c2" ,
                              "c1_name" , "c2_name" , "sigma","df",
                              "estimate","statistic","p.value",
                              "conf.low","conf.high","FDR")

# ContrastSimpleImpute----
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
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua::data_ionstar$normalized()
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
#'
#' pl <- tmp$get_Plotter()
#' pl$histogram()
#' pl$histogram_estimate()
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
    #' @description initialize
    #' @param lfqdata LFQData
    #' @param contrasts array of contrasts (see example)
    #' @param modelName default "groupAverage"
    initialize = function(lfqdata,
                          contrasts,
                          confint = 0.95,
                          p.adjust = prolfqua::adjust_p_values,
                          modelName = "groupAverage",
                          method = "V1",
                          probs = 0.1,
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
    #' @description get contrasts sides
    #'
    get_contrast_sides = function(){
      if (is.null(self$contrast_result)) {
        self$get_contrasts()
      }
      self$contrast_result %>% dplyr::select(contrast,c1 = c1_name,c2 = c2_name) %>% distinct()
    },
    #' @description table with results of contrast computation
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

        result <- result %>% rename(estimate = estimate_median, sigma = sd)
        result <- mutate(result, modelName = self$modelName, .before = 1)
        self$contrast_result <- result
      }
      invisible(ungroup(self$contrast_result))
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
        self$get_contrasts(),
        subject_Id = self$subject_Id,
        volcano = list(list(score = "FDR", thresh = 0.1)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "estimate",
        contrast = "contrast")
      return(res)
    },
    #' @description convert contrast results to wide format
    #' @param columns value column default p.value
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("estimate"),
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
#' library(prolfqua)
#' library(tidyverse)
#'
#' istar <- prolfqua::data_ionstar$normalized()
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
#'     "dil.e_vs_a" = "dilution.e - dilution.b" )
#'  #prolfqua::Contrasts$debug("get_linfct")
#'  #Contrasts$debug("get_contrasts")
#' contrastX <- prolfqua::Contrasts$new(mod,
#'  Contr)
#'
#' contrastX$get_contrast_sides()
#'
#' contrastX$get_linfct()
#' xx <- contrastX$get_contrasts()
#'
#'
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
    #' @description initialize
    #' create Contrast
    #' @param model a dataframe with a structure similar to that generated by \code{\link{build_model}}
    #' @param contrasts a character vector with contrast specificiation
    #' @param p.adjust function to adjust the p-values
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
    #' @description get both sides of contrasts
    get_contrast_sides = function(){
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
    #' @return data.frame with contrasts
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

        contrast_result <- self$p.adjust(contrast_result,
                                         column = "p.value",
                                         group_by_col = "contrast",
                                         newname = "FDR")
        contrast_result <- mutate(contrast_result, modelName = self$modelName, .before = 1)
        self$contrast_result <- dplyr::ungroup(contrast_result )
      }

      res <- if (!all) {
        self$contrast_result %>%
          select( -all_of(c("sigma.model",
                            "df.residual.model",
                            "std.error",
                            "isSingular")))

      }else{
        self$contrast_result
      }
      return(ungroup(res))
    },
    #' @description write results
    #' @param path directory
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write = function(path, filename, format = "xlsx"){
      filename <- if (missing(filename)) {self$modelName} else (filename )

      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",filename),
                      format = format)
    },
    #' @description
    #' return Contrast_Plotter
    get_Plotter = function(){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        volcano = list(list(score = "FDR", thresh = 0.1)),
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
#' library(prolfqua)
#' istar <- prolfqua::data_ionstar$normalized()
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
#'  contrast <- ContrastsModerated$new(contrast)
#'  bb <- contrast$get_contrasts()
#'  bb
#'  plotter <- contrast$get_Plotter()
#'  plotter$histogram()
#'  plotter$ma_plot()
#'  plotter$volcano()
#'
#' #bb %>% dplyr::rename(log2FC = estimate, mean_c1 = c1, mean_c2 = c2)
ContrastsModerated <- R6::R6Class(
  classname = "ContrastsModerated",
  inherit = ContrastsInterface,
  public = list(
    Contrast = NULL,
    #' @description initialize
    #' @param Contrast class implementing a method get_contrasts e.g. Contrasts
    modelName = character(),
    subject_Id = character(),
    p.adjust = NULL,
    initialize = function(Contrast,
                          modelName = paste0(Contrast$modelName, "_moderated"),
                          p.adjust = prolfqua::adjust_p_values
    ){
      self$Contrast = Contrast
      self$subject_Id = Contrast$subject_Id
      self$modelName = modelName
      self$p.adjust = p.adjust
    },
    #' @description get both sides of contrasts
    get_contrast_sides = function(){
      self$Contrast$get_contrast_sides()
    },
    #' @description get linear functions from contrasts
    #' @param global logical TRUE - get the a linear functions for all models, FALSE - linear function for each model
    get_linfct = function(global = TRUE){
      self$Contrast$get_linfct()
    },
    #' @description applies limma moderation
    #' @seealso \code{\link{moderated_p_limma_long}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      contrast_result <- self$Contrast$get_contrasts(all = FALSE)
      contrast_result <- moderated_p_limma_long(contrast_result ,
                                                group_by_col = "contrast")
      if (!all) {
        contrast_result <- contrast_result %>% select(-c( "sigma","df",
                                                          "statistic", "p.value","conf.low","conf.high",
                                                          "FDR",  "moderated.df.prior" ,
                                                          "moderated.var.prior"))
        contrast_result <- contrast_result %>% mutate(sigma = sqrt(moderated.var.post),.keep = "unused")
        contrast_result <- contrast_result %>% rename(
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
      return(dplyr::ungroup(contrast_result))
    },
    #' @description get \code{\link{Contrast_Plotter}}
    get_Plotter = function(){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        volcano = list(list(score = "FDR", thresh = 0.1)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "estimate",
        contrast = "contrast")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default moderated.p.value
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("estimate", columns),
                                                     contrast = 'contrast')
    },
    #' @description write results
    #' @param path directory
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write = function(path, filename, format = "xlsx"){
      filename <- if (missing(filename)) {self$modelName} else (filename )
      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",filename),
                      format = format)
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
#' istar <- prolfqua::data_ionstar$normalized()
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
#'  Contr <- c("dil.b_vs_a" = "dilution.a - dilution.b")
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
    Contrast = NULL,
    contrast_result = NULL,
    modelName = character(),
    subject_Id = character(),
    p.adjust = NULL,
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
    get_contrast_sides = function(){
      if (is.null(self$contrast_result)) {
        self$get_contrasts()
      }
      self$contrast_result %>% dplyr::select(contrast,c1 = c1_name,c2 = c2_name) %>% distinct()
    },
    get_linfct = function(){
      self$contrast$get_linfct()
    },

    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      if (is.null(self$contrast_result)) {
        contrast_result <- self$Contrast$get_contrasts(all = FALSE)
        contrast_result <- summary_ROPECA_median_p.scaled(
          contrast_result,
          contrast = "contrast",
          subject_Id = self$subject_Id[length(self$subject_Id) - 1],
          estimate = "estimate",
          statistic = "statistic",
          p.value = "p.value",
          max.n = 10)
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
    #' @description get \code{\link{Contrast_Plotter}}
    get_Plotter = function(){
      contrast_result <- self$get_contrasts()
      res <- Contrasts_Plotter$new(
        contrast_result,
        subject_Id = self$subject_Id[1],
        volcano = list(list(score = "FDR.beta.based.significance", thresh = 0.1)),
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
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(
        contrast_minimal,
        subject_Id = self$subject_Id[length(self$subject_Id) - 1],
        columns = c("estimate", columns),
        contrast = 'contrast')
      return(contrasts_wide)
    },
    #' @description write results
    #' @param path directory
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write = function(path, filename, format = "xlsx"){
      filename <- if (missing(filename)) {self$modelName} else (filename )

      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",filename),
                      format = format)
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
    contrast_result = NULL,
    modelName = character(),
    subject_Id = character(),
    initialize = function(contrastsdf,
                          subject_Id = "Prey",
                          modelName = "ContrastSaint"
    ){


      self$contrast_result = contrastsdf
      self$subject_Id = subject_Id
      self$modelName = modelName

      self$contrast_result <- contrastsdf %>% mutate(log2FC = log2(FoldChange),
                                                     c1_name = "Control",
                                                     c2_name = Bait,
                                                     c1 = log2(AvgIntensity) - log2(FoldChange)/2,
                                                     c2 = log2(AvgIntensity) + log2(FoldChange)/2 ,
                                                     modelName = modelName)

    },
    get_contrast_sides = function(){
      self$contrast_result %>% dplyr::select(Bait,c1 = c1_name,c2 = c2_name) %>% distinct()
    },
    get_linfct = function(){
      NULL
    },

    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)
    get_contrasts = function(all = FALSE){
      res <- self$contrast_result %>% select(all_of(c(self$subject_Id,
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
    #' @description get \code{\link{Contrast_Plotter}}
    #' @param fcthreshold fold change threshold to show
    #' @param scthreshold BFDR threshold to show in the heatmap.
    get_Plotter = function(fcthreshold = 1, bfdrthreshold = 0.1){
      res <- Contrasts_Plotter$new(
        self$contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = fcthreshold,
        volcano = list(list(score = "BFDR", thresh = bfdrthreshold)),
        histogram = list(list(score = "BFDR", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "log2FC",
        contrast = "Bait")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default SaintScore, BFDR
    to_wide = function(columns = c("SaintScore", "BFDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("log2FC", columns),
                                                     contrast = 'Bait')
      return(contrasts_wide)
    },
    #' @description write results
    #' @param path directory
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write = function(path, filename, format = "xlsx"){
      filename <- if (missing(filename)) {self$modelName} else (filename )
      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",filename),
                      format = format)
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
#' library(prolfqua)
#' library(tidyverse)
#' bb <- prolfqua::data_ionstar$normalized()
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
    contrast_result = NULL,
    modelName = character(),
    subject_Id = character(),
    initialize = function(contrastsdf,
                          subject_Id = "protein_Id",
                          modelName = "ContrastTable"
    ){
      self$contrast_result = contrastsdf
      self$subject_Id = subject_Id
      self$modelName = modelName
    },
    get_contrast_sides = function(){
      self$contrast_result %>% dplyr::select(contrast,c1 = c1_name,c2 = c2_name) %>% distinct()
    },
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
    #' @description get \code{\link{Contrast_Plotter}}
    get_Plotter = function(fcthreshold = 1, fdrthreshold = 0.1){
      res <- Contrasts_Plotter$new(
        self$contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = fcthreshold,
        volcano = list(list(score = "FDR", thresh = fdrthreshold)),
        histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                         list(score = "FDR", xlim = c(0,1,0.05))),
        modelName = self$modelName,
        estimate = "estimate",
        contrast = "contrast")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default beta.based.significance
    to_wide = function(columns = c("p.value", "FDR")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(contrast_minimal,
                                                     subject_Id = self$subject_Id,
                                                     columns = c("estimate", columns),
                                                     contrast = 'contrast')
      return(contrasts_wide)
    },
    #' @description write results
    #' @param path directory
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write = function(path, filename, format = "xlsx"){
      filename <- if (missing(filename)) {self$modelName} else (filename )
      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",filename),
                      format = format)
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
#' library(prolfqua)
#' library(tidyverse)
#'
#' istar <- prolfqua::data_ionstar$normalized()
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
#'   "dil.e_vs_a" = "dilution.e - dilution.a",
#'   "dil.e_vs_b" = "dilution.e - dilution.b",
#'   "dil.c_vs_b" = "dilution.c - dilution.b"
#'  )
#' #Contrasts$debug("get_contrasts")
#' contrast <- prolfqua::Contrasts$new(mod,
#'   Contr)
#' tmp <- contrast$get_contrasts()
#'
#' cp <- Contrasts_Plotter$new(tmp ,
#'  contrast$subject_Id,
#' volcano = list(list(score = "FDR", thresh = 1)),
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
    fcthresh = 1,
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
                          volcano = list(list(score = "FDR", thresh = 0.1)),
                          histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
                                           list(score = "FDR", xlim = c(0,1,0.05))),
                          fcthresh = 1,
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
      self$fcthresh = fcthresh
      self$contrast = contrast
    },
    #' @description  plot histogram of selected scores (e.g. p-value, FDR, t-statistics)
    histogram = function(){
      res <- list()
      if (length(self$histogram_spec) > 0) {
        for (spec in self$histogram_spec) {
          fig <- private$.histogram(score = spec )
          res[[spec$score]] <- fig
          self$figures[[paste0("histogram_", spec$score)]] <-
            list(fig = fig,
                 name = paste0(self$prefix, "_Histogram_", spec$score, "_", self$modelName))
        }
        return(res)
      }
    },
    #' @description plot histogram of estimated fold change
    #' @param binwidth with of bin in histogram
    histogram_estimate = function(binwidth = 0.05){
      re <- range(self$contrastDF[[self$estimate]], na.rm = TRUE)
      re[1] <- floor(re[1])
      re[2] <- ceiling(re[2])
      fig <- private$.histogram(score = list(score =  self$estimate, xlim = c(re,binwidth)))
      self$figures[["histogram_estimate"]] <- list(fig = fig,
                                                   name = paste0(self$prefix,"_Histogram_Estimate_", self$modelName))
      return(fig)

    },
    #' @description volcano plots (fold change vs FDR)
    volcano = function(colour = "modelName"){
      fig <- private$.volcano(self$contrastDF, self$volcano_spec, colour = colour )
      self$figures[["volcano"]] <- list(fig = fig, name = paste0(self$prefix, "_Volcano_", self$modelName))
      return(fig)
    },
    #' @description plotly volcano plots
    volcano_plotly = function(colour = "modelName"){
      contrastDF <- self$contrastDF %>% plotly::highlight_key(~ subject_Id)
      res <- private$.volcano(contrastDF, self$volcano_spec, colour = colour )
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
    ma_plot = function(fc = 1, colour = "modelName"){
      contrastDF <- self$contrastDF
      if (!is.null(contrastDF$c1) && !is.null(contrastDF$c2)) {
        # pdf version
        fig <- private$.ma_plot(contrastDF ,self$contrast, colour = colour)
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
    ma_plotly = function(fc =1, colour = "modelName"){
      # html version
      contrastDF  <- self$contrastDF
      if (!is.null(contrastDF$c1) && !is.null(contrastDF$c2)) {

        contrastDF  <- contrastDF %>%
          plotly::highlight_key(~subject_Id)
        fig_plotly <- private$.ma_plot(contrastDF, self$contrast, colour = colour) %>%
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
    .volcano = function(contrasts, scores,  colour = NULL){
      fig <- list()
      for (score in scores) {
        column <- score$score
        fig[[column]] <- prolfqua:::.multigroupVolcano(
          contrasts,
          effect = self$estimate,
          p.value = column,
          condition = self$contrast,
          text = "subject_Id",
          xintercept = c(-self$fcthresh, self$fcthresh),
          yintercept = score$thresh,
          colour = colour,
          scales = "free_y")

      }
      return(fig)

    },
    .histogram  = function(score){
      xlim = score$xlim
      score = score$score
      plot <- self$contrastDF %>% ggplot(aes(x = !!sym(score))) +
        geom_histogram(breaks = seq(from = xlim[1], to = xlim[2], by = xlim[3])) +
        facet_wrap(vars(!!sym(self$contrast)))
      return(plot)
    },
    .ma_plot = function(x, contrast, colour = NULL){
      x <- ggplot(x , aes(x = (c1 + c2)/2,
                          y = !!sym(self$estimate),
                          text = !!sym("subject_Id"),
                          colour = !!sym(colour))) +
        geom_point(alpha = 0.5) +
        scale_colour_manual(values = c("black", "green")) +
        facet_wrap(vars(!!sym(contrast))) + theme_light() +
        geom_hline(yintercept = c(-self$fcthresh, self$fcthresh), linetype = "dashed",colour = "red") +
        ylab("log fold change (M)") + xlab("mean log intensities (A)")
      return(x)
    }
  )
)




# Merge contrasts ----
#' add contrast results from two different functions. Tipically used with Contrast and Cotnrast simple imputed.
#'
#' @param prefer contrasts to use preferentially
#' @param add contrasts to add from if missing in prefer
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


