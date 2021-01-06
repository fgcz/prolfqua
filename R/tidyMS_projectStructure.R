.checkForFile <- function(inputData){
  if (!is.null(inputData)) {
    if ( file.exists(inputData) ) {
      inputData
    } else {
      stop("File does not exist : ", inputData)
    }
  }
}

.dirmaker <- function(path){
  if (!dir.exists(path)) {
    dir.create(path)
  }
  NULL
}

.dircleaner <- function(path){
  if (unlink(path , recursive = TRUE) != 0 ) {
    message("could not clean : ", path)
  }
  NULL
}

#' keep track of folder paths and create them if needed
#'
#' @export
#'
#' @examples
#' tmp <- ProjectStructure$new("./test_project",
#' project_Id  = 3000,
#' order_Id = 6200,
#' workunit_Id = 23000,
#' inputAnnotation = ".",
#' inputData = "."
#' )
#' tmp$qc_path()
#' tmp$modelling_path()
#'
#' tmp$modelling_path()
#' tmp$modelling_dir
#' tmp$modelling_path("second_model")
#' tmp$create()
#' tmp$reset()
#'
ProjectStructure <-
  R6::R6Class("ProjectStructure", public = list(
    #' @field outpath path
    #' @field qc_dir qc results directory name
    #' @field modelling_dir modeling results directory name
    #' @field project_Id project_Id
    #' @field order_Id order_Id
    #' @field workunit_Id workunit_Id
    #' @field inputData inputFile
    #' @field inputAnnotation inputAnnotation xlsx
    outpath = "",
    qc_dir = "",
    modelling_dir = "",
    project_Id = integer(),
    order_Id = integer(),
    workunit_Id = integer(),
    inputData = character(),
    inputAnnotation = character(),
    #' @description
    #' create ProjectStructure
    #' @param outpath directory
    #' @param project_Id bfabric project ID
    #' @param workunit_Id bfabric workunit_Id
    #' @param order_Id bfabric order_Id
    #' @param inputAnnotation input annotation path
    #' @param inputData input data path
    #' @param qc_dir qc folder
    #' @param modelling_dir modelling results folder

    initialize = function(outpath,
                          project_Id,
                          order_Id,
                          workunit_Id,
                          inputAnnotation,
                          inputData,
                          qc_dir = "qc_results",
                          modelling_dir = "modelling_results"){
      self$outpath = outpath
      self$project_Id = project_Id
      self$order_Id = order_Id
      self$workunit_Id = workunit_Id
      self$inputData = .checkForFile(inputData)
      self$inputAnnotation = .checkForFile(inputAnnotation)
      self$qc_dir = qc_dir
      self$modelling_dir = modelling_dir
    },
    #' @description
    #' create outpath
    create_outpath = function(){
      .dirmaker(self$outpath)
    },
    #' @description
    #' create qc dir
    #' @param qc_dir QC directory
    qc_path = function(qc_dir){
      if (!missing(qc_dir)) {
        self$qc_dir <- c(self$qc_dir, qc_dir)
      }
      qcpath <- file.path(self$outpath, self$qc_dir)
      qcpath
    },
    #' @description
    #' create modelling path
    #' @param modelling_dir directory with modelling data
    modelling_path = function(modelling_dir){
      if (!missing(modelling_dir)) {
        self$modelling_dir <- c(self$modelling_dir, modelling_dir)
      }
      modellingpath <- file.path(self$outpath, self$modelling_dir)
      modellingpath
    },
    #' @description
    #' create all directories
    create = function(){
      .dirmaker(self$outpath)
      sapply(self$qc_path(),.dirmaker)
      sapply(self$modelling_path(),.dirmaker)
      NULL
    },
    #' @description
    #' empty modelling_path and qc_path folder.
    reset = function(){
      sapply(self$qc_path(), .dircleaner)
      sapply(self$modelling_path(), .dircleaner)
      NULL
    }

  )
  )

