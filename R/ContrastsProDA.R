#' ContrastsProDA Wrapper to results produced by proDA
#'
#' @export
#' @family modelling
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config()
#' istar$config <- istar$config
#' istar_data <- istar$data
#' lfd <- LFQData$new(istar_data, istar$config)
#' se <- prolfqua::LFQDataToSummarizedExperiment(lfd)
#' if(require(proDA)){
#' fit <- proDA::proDA(se, design = ~ group_ - 1, data_is_log_transformed = TRUE)
#' contr <- list()
#' contrasts <- c("group_AvsCtrl" = "group_A - group_Ctrl",
#'                "group_BvsCtrl" = "group_B - group_Ctrl")
#' contr[["group_AvsCtrl"]] <- data.frame(
#'   contrast = "group_AvsCtrl",
#'   proDA::test_diff(fit, contrast = "group_A - group_Ctrl"))
#' contr[["group_BvsCtrl"]] <- data.frame(
#'   contrast = "group_BvsCtrl",
#'   proDA::test_diff(fit, contrast = "group_B - group_Ctrl"))
#'
#' bb <- dplyr::bind_rows(contr)
#' cproDA <- ContrastsProDA$new(bb, contrasts = contrasts, subject_Id = "name")
#' x <- cproDA$get_contrasts()
#' cproDA$get_linfct()
#' contsides <- cproDA$get_contrast_sides()
#' stopifnot(ncol(cproDA$to_wide()) == c(7))
#' tmp <- cproDA$get_Plotter()
#' tmp$volcano()$pval
#' tmp$volcano()$adj_pval
#' }
#'
ContrastsProDA <- R6::R6Class(
  "ContrastsProDA",
  inherit = prolfqua::ContrastsInterface,
  public = list(
    #' @field  contrast_result contrast result
    contrast_result = NULL,
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id columns with protein ID's
    subject_Id = character(),
    #' @field contrasts named vector of length 1.
    contrasts = character(),
    #' @description
    #' initialize
    #' @param contrastsdf data.frame returned by proDA
    #' @param contrasts  contrasts
    #' @param subject_Id column name with protein ID's
    #' @param modelName name of model default value ContrastProDA
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
    #' @description
    #' show names of contrasts
    #' @return data.frame
    get_contrast_sides = function(){
      # extract contrast sides
      tt <- self$contrasts[grep("-",self$contrasts)]
      tt <- tibble(contrast = names(tt) , rhs = tt)
      tt <- tt |> mutate(rhs = gsub("[` ]","",rhs)) |>
        tidyr::separate(rhs, c("group_1", "group_2"), sep = "-")
      return(tt)
    },
    #' @description
    #' get linear function used to determine contrasts
    #' @return data.frame
    get_linfct = function(){
      NULL
    },
    #' @description
    #' get contrasts
    #' @param all (default FALSE)
    get_contrasts = function(all = FALSE){
      return(self$contrast_result)
    },
    #' @description get \code{\link{Contrast_Plotter}}
    #' @param fcthreshold fold change threshold to show
    #' @param fdrthreshold FDR threshold
    #' @param tstatthreshold t statistics threshold
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
    #' @param columns value column default t_statistic, adj_pval
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
    #' @param filename file to write to
    #' @param format default xlsx \code{\link{lfq_write_table}}
    write = function(path, filename, format = "xlsx"){
      filename <- if (missing(filename)) {self$modelName} else (filename )
      lfq_write_table(self$get_contrasts(),
                      path = path,
                      name  = paste0("Contrasts_",filename),
                      format = format)
    }
  ))
