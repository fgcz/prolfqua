#' find file stored in package
#' @export
#' @keywords internal
#' @examples
#' find_package_file("prolfqua","extdata/medata.csv")
#'
find_package_file <- function(packagename, file){
  src_script <- file.path(find.package(packagename) , file )
  if (!file.exists(src_script)) {
    src_script <- file.path(find.package(packagename) , "inst" , file)
  }
  if (file.exists(src_script)) {
    return(src_script)
  }else{
    return(NULL)
  }
}

#' copy script files and other from a package to workdir
#' @export
#' @keywords internal
scriptCopyHelperVec <-
  function(runscripts,
           workdir = getwd(),
           packagename = "prolfqua") {
    res <- NULL
    for (scripts in runscripts) {
      src_script <- file.path(find.package(packagename) , scripts)
      dest_script <- file.path(workdir , basename(scripts))
      message("copy ", src_script, " to ", dest_script)
      if (!file.exists(src_script)) {
        src_script <- file.path(find.package(packagename) , "inst" , scripts)
        if (!file.exists(src_script)) {
          warning(paste("could not copy script file.", dest_script, sep = " "))
        }
      }
      if (!file.copy(src_script , dest_script, overwrite = TRUE)) {
        warning(paste("could not copy script file.", src_script, " to ", dest_script, sep = " "))
      } else{
        res <- c(res, dest_script)
      }
    }
    message(paste(
      "your working directory now should contain: ",
      length(res) ,
      "new files :\n",
      sep = " "
    ))
    return(res)
  }


.run_markdown_with_params <-
  function(params,
           markdown_path,
           dest_path,
           dest_file_name,
           workdir = tempdir(),
           packagename = "prolfqua",
           format = "pdf") {
    res <- prolfqua::scriptCopyHelperVec(markdown_path,
                                workdir = workdir,
                                packagename = packagename)
    dist_file_path <-
      file.path(dest_path, paste0(dest_file_name, ".", format))
    if (is.null(res)) {
      return(NULL)
    }
    rmarkdown::render(
      res[1],
      output_format = if (format == "pdf") {
        bookdown::pdf_document2()
      } else{
        bookdown::html_document2()
      },
      params = params,
      envir = new.env()
    )

    pdf_doc <- paste0(tools::file_path_sans_ext(res[1]), ".", format)
    message("XXXX--------------------------------------XXXX")
    if (pdf_doc != dist_file_path) {
      message("from " , pdf_doc, " to ", dist_file_path)
      file.copy(pdf_doc, dist_file_path, overwrite = TRUE)
    }
    return(dist_file_path)
  }

#' render MQ Summary.
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param project_config list with workunit_Id project_Id order_Id
#' @param pep are these peptide or protein data
#' @param dest_path destination path
#' @param dest_file_name name of pdf file
#' @param workdir working directory
#' @param format either pdf or html
#' @family vignetteHelpers
#' @keywords internal
#' @export
#' @examples
#' bb <- prolfqua_data('data_skylinePRMSample_A')
#' config <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' projectConfig <- list(workunit_Id = "xx", project_Id = "xy", order_Id = "z")
#' if(FALSE){
#' render_MQSummary_rmd(analysis,
#'   config ,
#'   projectConfig,
#'   workdir= tempdir()) # tempdir(check = FALSE))
#' }
render_MQSummary_rmd <-
  function(pdata,
           config,
           project_conf,
           pep = TRUE,
           dest_path = ".",
           dest_file_name = "QCandSampleSize.Rmd",
           workdir = tempdir(),
           format = c("pdf", "html"),
           markdown_path = c("doc/QCandSampleSize.Rmd"))
  {
    dist_file_path <- .run_markdown_with_params(
      list(
        data = pdata,
        configuration = config,
        project_conf = project_conf,
        pep = pep
      ),
      markdown_path = markdown_path,
      dest_path = dest_path,
      dest_file_name = dest_file_name,
      workdir = workdir,
      packagename = "prolfqua",
      format = match.arg(format)
    )
    return(dist_file_path)
  }



