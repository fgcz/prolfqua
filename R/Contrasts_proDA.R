# ContrastsSaintExpress ----
#' ContrastsSaintExpress Wrapper to results produced by SaintExpress (list.txt file)
#'
#' @export
#' @family modelling
#'
Contrasts_proDA <- R6::R6Class(
  "Contrasts_proDA",
  inherit = prolfqua::ContrastsInterface,
  public = list (
    contrast_result = NULL,
    modelName = character(),
    subject_Id = character(),
    #' @field contrast named vector of length 1.
    contrast = character(),
    initialize = function(contrastsdf,
                          contrast,
                          subject_Id = "protein_Id",
                          modelName = "ContrastSaint"
    ){

      self$contrast_result = contrastsdf
      self$subject_Id = subject_Id
      self$modelName = modelName
      self$contrast = contrast


      contrastdf$contrast <- name( contrast )
      cs <- .get_sides(contrast)
      contrastdf$g1_name <- cs[1]
      contrastdf$g2_name <- cs[2]


      self$contrast_result <- contrastsdf

      mapping <- c(subject_Id = "name",
        p.value = "pval",
        FDR = "adj_pval",
        diff  = "diff",
        statistic = "t_statistic",
        se = "se",
        df = "df",
        avg_abundance = "avg_abundance",
        n_approx = "n_approx",
        n_obs = "n_obs")
    },
    get_contrast_sides = function(){
      res <- data.frame(contrast = contrast, g1 = cs[1], g2 = cs[2])
      return(res)
    },
    get_linfct = function(){
      # maybe later.
      NULL
    },

    add = function(...){
      all <- list(...)
      all[[length(all)+1]] <- self

    },
    #' @description
    #' get contrasts
    #' @seealso \code{\link{summary_ROPECA_median_p.scaled}}
    #' @param all should all columns be returned (default FALSE)
    #' @param global use a global linear function (determined by get_linfct)

    get_contrasts = function(all = FALSE){
      res <- self$contrast_result
      res
    },
    #' @description get \code{\link{Contrast_Plotter}}
    #' @param fcthreshold fold change threshold to show
    #' @param scthreshold BFDR threshold to show in the heatmap.
    get_Plotter = function(fcthreshold = 1, bfdrthreshold = 0.1){
      res <- ContrastsPlotter$new(
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
