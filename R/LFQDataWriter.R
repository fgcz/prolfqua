# LFQDataWriter -----
#'
#' Write LFQData, or provide outputs for writing.
#'
#' @family LFQData
#'
#' @seealso \code{\link{LFQData}}
#' @export
#'
LFQDataWriter <- R6::R6Class(
  "LFQDataWriter",list(
    #' @field lfq LFQData
    #' @field format format to write to
    #' @field prefix prefix of filename
    #' @field file_paths list with paths were data was written to.
    #'
    lfq = NULL,
    format = "",
    prefix = "",
    file_paths = list(),
    #' @description
    #' initialize class
    #' @param lfqdata LFQData
    #' @param prefix prefix files with
    #' @param format which format to write to ("xlsx", "csv", "html")
    initialize = function(lfqdata,  prefix = "ms_", format="xlsx"){
      self$lfq = lfqdata
      self$format = format
      self$prefix = prefix
    },
    #' @description
    #' Get data in long format for writing
    #' @return tibble
    get_long = function(){
      #' gets data formatted for writing
      separate_factors(
        separate_hierarchy(self$lfq$data, self$lfq$config),
        self$lfq$config)
    },
    #' @description
    #' Get data in Wide format for writing
    #' @return list with data and annotation
    get_wide = function(){
      #' gets data formatted for writing
      wide <- self$lfq$to_wide()
      res <- list(data = separate_hierarchy(wide$data, self$lfq$config), annotation = wide$annotation)
      return(res)
    },
    #' @description
    #' write data to file
    #' @param path_qc path to write to
    write_long = function(path_qc) {
      fname <- paste0(self$prefix,"intensities_long")
      self$file_paths[[fname]] <-
        lfq_write_table(self$get_long(),
                        path = path_qc,
                        name = fname,
                        format = self$format)
    },
    #' @description
    #' write data to file
    #' @param path_qc path to write to
    write_wide = function(path_qc) {

      wide <- self$get_wide()
      fname <- paste0(self$prefix,"intensities_wide")
      self$file_paths[[fname]] <-
        lfq_write_table(wide$data,
                        path = path_qc,
                        name = fname,
                        format = self$format)

      fname <- paste0(self$prefix,"intensities_file_annotation")
      self$file_paths[[fname]] <-
        lfq_write_table(wide$annotation,
                        path = path_qc,
                        name = fname,
                        format = self$format)
    }
  ))
