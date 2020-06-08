ProjectStructure <-
  R6::R6Class("ProjectStructure", public = list(
    outpath = "",
    qc_path = "",
    modelling_path = "",
    #' @description
    #' create ProjectStructure
    #' @param outpath directory
    #' @param qc_path qc folder
    #' @param modelling_path modelling results folder
    initialize = function(outpath ,
                          qc_path = "qc_results",
                          modelling_path = "modelling_results"){
      self$outpath = outpath
      self$qc_path = file.path(outpath, qc_path )
      self$modelling_path = file.path(outpath, modelling_path )
    },
    create_outpath = function(){
      if (!dir.exists(self$outpath)) {
        dir.create(self$outpath)
      }
    },
    creat_qc_dir = function(){
      self$create_outpath()
      if (!dir.exists(self$qc_path)) {
        dir.create(self$qc_path)
      }
    },
    create_modelling_path = function(){
      self$create_outpath()
      if (!dir.exists(self$modelling_path)) {
        dir.create(self$modelling_path)
      }
    },
    create = function(){
      self$creat_qc_dir()
      self$create_modelling_path()
    }

  )
  )

#' Create Project structure
#' @export
#' @examples
#' tmp <- createProjectStructure("./test_project")
#' tmp$qc_path
#' tmp$modelling_path
#' if(FALSE){tmp$create()}
#'
createProjectStructure <- function(outpath,
                                   qc_path = "qc_results",
                                   modelling_path = "modelling_results"){
  res <- ProjectStructure$new(outpath, qc_path, modelling_path)
  return(res)
}
