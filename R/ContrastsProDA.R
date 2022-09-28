#' ContrastsProDA Wrapper to results produced by proDA
#'
#' @export
#' @family modelling
#' @examples
#'
#' istar <- old2new( prolfqua_data('data_ionstar')$normalized() )
#' istar_data <- dplyr::filter(istar$data ,protein_Id %in% sample(protein_Id, 10))
#' lfd <- LFQData$new(istar_data, istar$config)
#' se <- prolfqua::LFQDataToSummarizedExperiment(lfd)
#'
#' fit <- proDA::proDA(se, design = ~ dilution. - 1,data_is_log_transformed = TRUE)
#' contr <- list()
#'
#' contrasts <- c("dilution_(9/7.5)_1.2" = "dilution.e - dilution.d",
#'                "dilution_(7.5/6)_1.25" = "dilution.d - dilution.c")
#' contr[["dilution_(9/7.5)_1.2"]] <- data.frame(
#'   contrast = "dilution_(9/7.5)_1.2",
#'   proDA::test_diff(fit, contrast = "dilution.e - dilution.d"))
#' contr[["dilution_(7.5/6)_1.25"]] <- data.frame(
#'   contrast = "dilution_(7.5/6)_1.25",
#'   proDA::test_diff(fit, contrast = "dilution.d - dilution.c"))
#'
#' bb <- dplyr::bind_rows(contr)
#' cproDA <- ContrastsProDA$new(bb, contrasts = contrasts, subject_Id = "name")
#' cproDA$get_contrasts() |> head()
#' cproDA$get_linfct()
#' contsides <- cproDA$get_contrast_sides()
#' stopifnot(ncol(cproDA$to_wide()) == c(7))
#' tmp <- cproDA$get_Plotter()
#' tmp$volcano()$pval
#' tmp$volcano()$adj_pval

ContrastsProDA <- R6::R6Class(
  "ContrastsProDA",
  inherit = prolfqua::ContrastsInterface,
  public = list (
    contrast_result = NULL,
    modelName = character(),
    subject_Id = character(),
    #' @field contrast named vector of length 1.
    contrasts = character(),
    initialize = function(contrastsdf,
                          contrasts,
                          subject_Id = "name",
                          modelName = "ContrastProDA"
    ){
      if(!"modelName" %in% names(contrastsdf)){
        contrastsdf$modelName <- modelName
      }
      self$contrast_result = contrastsdf
      self$subject_Id = subject_Id
      self$modelName = modelName
      self$contrasts = contrasts
      self$contrast_result <- contrastsdf

    },
    get_contrast_sides = function(){
      # extract contrast sides
      tt <- self$contrasts[grep("-",self$contrasts)]
      tt <- tibble(contrast = names(tt) , rhs = tt)
      tt <- tt |> mutate(rhs = gsub("[` ]","",rhs)) |>
        tidyr::separate(rhs, c("group_1", "group_2"), sep = "-")
      return(tt)
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
      return(self$contrast_result)
    },
    #' @description get \code{\link{Contrast_Plotter}}
    #' @param fcthreshold fold change threshold to show
    #' @param scthreshold BFDR threshold to show in the heatmap.
    get_Plotter = function(fcthreshold = 1, fdrthreshold = 0.1, tstatthreshold = 5){
      res <- ContrastsPlotter$new(
        self$contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = fcthreshold,
        volcano = list(list(score = "pval", thresh = fdrthreshold),list(score = "adj_pval", thresh = fdrthreshold)),
        histogram = list(list(score = "pval", xlim = c(0,1,0.05)), list(score = "adj_pval", xlim = c(0,1,0.05))),
        score = list(list(score = "t_statistic", thresh = tstatthreshold )),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast")
      return(res)
    },
    #' @description convert to wide format
    #' @param columns value column default SaintScore, BFDR
    to_wide = function(columns = c("t_statistic", "adj_pval")){
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(
        contrast_minimal,
        subject_Id = self$subject_Id,
        columns = c("diff", columns),
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
