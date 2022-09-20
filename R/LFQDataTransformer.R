# LFQDataTransformer ----

#' Decorate LFQData with Methods for transforming Intensities
#'
#' @export
#'
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
#' lfqdata <- LFQData$new(data, istar$config)
#'
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#'
#' x <- lfqTrans$intensity_array(log2)
#' x$lfq$config$table$is_intensity_transformed
#' x <- x$intensity_matrix(robust_scale)
#' plotter <- x$lfq$get_Plotter()
#' plotter$intensity_distribution_density()
#'
#' # transform by asinh root and scale
#'
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' x <- lfqTrans$intensity_array(asinh)
#' mads1 <- mean(x$get_scales()$mads)
#' x <- lfqTrans$intensity_matrix(robust_scale, force = TRUE)
#' mads2 <- mean(x$get_scales()$mads)
#'
#' stopifnot(abs(mads1 - mads2) < 1e-8)
#'
#'
#' stopifnot(abs(mean(x$get_scales()$medians)) < 1e-8)
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' lfqTrans$log2()
#' before <- lfqTrans$get_scales()
#' lfqTrans$robscale()
#' after <- lfqTrans$get_scales()
#' stopifnot(abs(mean(before$medians) - mean(after$medians)) < 1e-8)
#' stopifnot(abs(mean(before$mads) - mean(after$mads)) < 1e-8)
#'
#' # normalize data using vsn
#' lfqcopy <- lfqdata$get_copy()
#' lfqTrans <- lfqcopy$get_Transformer()
#' lfqTransCheck <- lfqcopy$get_Transformer()
#'
#' lfqTransCheck$log2()
#' lfqTransCheck$get_scales()
#' lfqTransCheck$lfq$get_Plotter()$intensity_distribution_density()
#'
#' if(require("vsn")){
#'  res <- lfqTrans$intensity_matrix( .func = vsn::justvsn)
#'  res$lfq$get_Plotter()$intensity_distribution_density()
#'  res$get_scales()
#' }
#'
#'
LFQDataTransformer <- R6::R6Class(
  "LFQDataTransformer",
  public = list(
    #' @field lfq LFQData object
    lfq = NULL,
    #' @description
    #' initialize
    #' @param lfqdata
    #' LFQData object to transform
    initialize = function(lfqdata){
      self$lfq = lfqdata$clone(deep = TRUE)
    },
    #' @description
    #' log2 transform data
    #' @param force if FALSE, then data already log2 transformed will not be transformed a second time. TRUE force log transformation.
    #' @return LFQDataTransformer
    log2 = function(force = FALSE){
      if (self$lfq$is_transformed() == FALSE | force ) {
        self$lfq$data  <-  prolfqua::transform_work_intensity(self$lfq$data ,
                                                              self$lfq$config, log2)
        self$lfq$is_transformed(TRUE)
      } else {
        warning("data already transformed. If you still want to log2 tranform, set force = TRUE")
      }
      invisible(self)

    },
    #' @description
    #' get mean and variance and standard deviation in each sample
    #' @return list with means and mads
    get_scales = function()
    {
      get_robscales(self$lfq$data, self$lfq$config)
    },
    #' @description
    #' robust scale data
    #' @param colname new name of transformed column
    #' @param preserveMean should original mean value be preserved TRUE, if FALSE then center at zero
    #' @return LFQDataTransformer (self)
    robscale = function(preserveMean = TRUE, colname = "transformedIntensity"){
      res <- self$robscale_subset(self$lfq, preserveMean = preserveMean, colname = colname)
      invisible(res)
    },
    #' @description
    #' log2 transform and robust scale data based on subset
    #' @param lfqsubset LFQData subset to use for normalization
    #' @param preserveMean should original mean value be preserved TRUE, if FALSE then center at zero
    #' @param colname - how to name the transformed intensities, default transformedIntensity
    #' @return LFQDataTransformer (self)
    #'
    robscale_subset = function(lfqsubset,
                               preserveMean = TRUE,
                               colname = "transformedIntensity"){
      message("data is : ",self$lfq$is_transformed())
      if (self$lfq$is_transformed() != lfqsubset$is_transformed()) {
        warning("the subset must have the same config as self")
        invisible(NULL)
      }
      scales <- prolfqua::scale_with_subset(self$lfq$data,
                                            lfqsubset$data,
                                            self$lfq$config,
                                            preserveMean = preserveMean)
      self$lfq$data  <- scales$data
      if (!is.null(colname)) {
        self$lfq$data <- self$lfq$data |>
          dplyr::rename(!!colname := self$lfq$config$table$getWorkIntensity())
        self$lfq$config$table$popWorkIntensity()
        self$lfq$config$table$setWorkIntensity(colname)
      }
      invisible(self)

    },
    #' @description
    #' Transforms intensities
    #' @param .func transformation function working with arrays e.g. log2, log10, asinh etc.
    #' @param force transformation on already transformed data.
    #'
    #' @return LFQDataTransformer (self)
    #'
    intensity_array = function(.func = log2, force = FALSE) {
      if (!self$lfq$is_transformed() | force) {
        .call <- as.list( match.call() )
        r <- prolfqua::transform_work_intensity(
          self$lfq$data,
          self$lfq$config,
          .func = .func,
          .funcname = deparse(.call$.func))
        self$lfq$data <- r
        self$lfq$is_transformed(TRUE)

      } else {
        warning("data already transformed. If you still want to log2 tranform, set force = TRUE")

      }
      invisible(self)
    },
    #' @description
    #' pass a function which works with matrices, e.g., vsn::justvsn
    #' @param .func any function taking a matrix and returning a matrix (columns sample, rows feature e.g. base::scale) default robust_scale
    #' @param force transformation on data already transformed
    #' @return LFQDataTransformer (self)
    #'
    intensity_matrix = function(.func = robust_scale, force = FALSE){
      if (!self$lfq$is_transformed() | force) {
        .call <- as.list( match.call() )
        r <- prolfqua::applyToIntensityMatrix(
          self$lfq$data,
          self$lfq$config,
          .func = .func,
          .funcname = deparse(.call$.func))
        self$lfq$data <- r
        self$lfq$is_transformed(TRUE)
      } else {
        warning("data already transformed. If you still want to log2 tranform, set force = TRUE")

      }
      invisible(self)
    }
  )
)
