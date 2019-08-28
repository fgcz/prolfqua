.scriptCopyHelperVec <- function(runscripts, workdir = getwd(), packagename = "LFQService" ){
  res <- NULL
  for(scripts in runscripts){
    src_script <- file.path( find.package(packagename) , scripts )
    dest_script <- file.path(workdir ,basename(scripts))
    message("copy ", src_script, " to ", dest_script)

    if(!file.copy(src_script , dest_script, overwrite = TRUE)){
      warning(paste("could not copy script file.", dest_script, sep=" "))
    }else{
      res <- c(res,dest_script )
    }
  }
  message(paste("your working directory now should contain: ", length(res) , "new files :\n",sep=" "))
  return(res)
}

.run_markdown_with_params <- function(params,
                                      markdown_path,
                                      dest_path,
                                      dest_file_name,
                                      workdir = tempdir(),
                                      packagename = "LFQService"){
  res <- .scriptCopyHelperVec(markdown_path, workdir = workdir, packagename = packagename)
  dist_file_path <- file.path(dest_path, dest_file_name)
  if(is.null(res)){
    return(NULL)
  }
  rmarkdown::render(res[1],
                    output_format = bookdown::pdf_document2(),
                    params=params,
                    envir = new.env()
  )

  pdf_doc <- paste0(tools::file_path_sans_ext(res[1]),".pdf")
  message("XXXX--------------------------------------XXXX")
  if(pdf_doc != dist_file_path){
    message("from " , pdf_doc, " to ", dist_file_path)
    file.copy(pdf_doc, dist_file_path, overwrite = TRUE)
  }
  return(dist_file_path)
}

#' render MQ Summary.
#' @export
#' @examples
#' render_MQSummary_rmd(LFQService::sample_analysis,LFQService::skylineconfig, workdir=tempdir(check = FALSE))
#'
render_MQSummary_rmd <- function(data, config,
                                 workunit_id = "w1",
                                 project_id = "p1000",
                                 order_id = "o2000",
                                 pep = TRUE,
                                 dest_path =".",
                                 dest_file_name = "MQSummary2.pdf",
                                 workdir = tempdir())
{
  dist_file_path <- .run_markdown_with_params(
    list(data = data,
         configuration=config,
         project_id = project_id,
         workunit_id = workunit_id,
         order_id = order_id,
         pep=pep),
    markdown_path =c("rmarkdown/MQSummary2.Rmd", "rmarkdown/CVReport.Rmd"),
    dest_path = dest_path,
    dest_file_name = dest_file_name,
    workdir = workdir,
    packagename = "LFQService"
  )
  return(dist_file_path)
}

#' render MQ Summary.
#' @export
#' @examples
#' #render_MQSummary_rmd(LFQService::sample_analysis,LFQService::skylineconfig, workdir=tempdir(check = FALSE))
#'
render_METABO_Summary_rmd <- function(data, config,
                                      dest_path =".",
                                      dest_file_name = "Metabo_Summary.pdf",
                                      workdir = tempdir())
{
  dist_file_path <- .run_markdown_with_params(
    list(data = data, configuration=config$clone(deep=TRUE)),
    markdown_path ="rmarkdown/METABO_Summary.Rmd",
    dest_path = dest_path,
    dest_file_name = dest_file_name,
    workdir = workdir,
    packagename = "LFQService"
  )
  return(dist_file_path)
}


#' render Filtering Summary.
#' @export
#' @examples
#'
render_SummarizeFiltering_rmd <- function(results,
                                          dest_path = ".",
                                          dest_file_name = "Summarize_Filtering.pdf",
                                          workdir = tempdir())
{
  dist_file_path <- .run_markdown_with_params(
    results,
    markdown_path ="rmarkdown/Summarize_Filtering.Rmd",
    dest_path = dest_path,
    dest_file_name = dest_file_name,
    workdir = workdir,
    packagename = "LFQService"
  )
}

#' render Filtering Summary.
#' @export
#' @examples
#'
render_METABO_SummarizeFiltering_rmd <- function(results,
                                                 dest_path = ".",
                                                 dest_file_name = "METABO_Summarize_Filtering.pdf",
                                                 workdir = tempdir())
{
  dist_file_path <- .run_markdown_with_params(
    results,
    markdown_path ="rmarkdown/METABO_Summarize_Filtering.Rmd",
    dest_path = dest_path,
    dest_file_name = dest_file_name,
    workdir = workdir,
    packagename = "LFQService"
  )
}

