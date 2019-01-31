.scriptCopyHelperVec <- function(runscripts, workdir = getwd() ){
  res <- NULL
  for(scripts in runscripts){
    src_script <- file.path( find.package("LFQService") , scripts )
    dest_script <- file.path(workdir ,basename(scripts))
    message("copy ", src_script, " to ", dest_script)

    if(!file.copy(src_script , dest_script, overwrite = TRUE)){
      warning(paste("could not copy script file.", dest_script, sep=" "))
    }else{
      res <- c(res,dest_script )
    }
  }
  message(paste("your working directory now should contain ", length(runscripts) , "new files :\n",sep=" "))
  return(res)
}

.run_markdown_with_params <- function(params,
                                     markdown_path,
                                     dest_path,
                                     dest_file_name,
                                     workdir = tempdir()){
  markdown_file <- basename(markdown_path)
  res <- .scriptCopyHelperVec(markdown_path, workdir = workdir)
  dist_file_path <- file.path(dest_path, dest_file_name)
  rmarkdown::render(res,
                    output_format = rmarkdown::pdf_document(),
                    params=params, envir = new.env()
  )

  pdf_doc <- paste0(tools::file_path_sans_ext(res),".pdf")
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
                              dest_path =".",
                              dest_file_name = "MQSummary.pdf",
                              workdir = tempdir())
  {
  dist_file_path <- .run_markdown_with_params(list(data = data, configuration=config$clone(deep=TRUE)),
                           markdown_path ="rmarkdown/MQSummary.Rmd",
                           dest_path = dest_path,
                           dest_file_name = dest_file_name,
                           workdir = workdir
                           )
  return(dist_file_path)
}

#' render MQ Summary.
#' @export
#' @examples
#'
render_SummarizeFiltering_rmd <- function(results,
                                          dest_path = ".",
                                          dest_file_name = "Summarize_Filtering.pdf",
                                          workdir = tempdir())
{
  dist_file_path <- .run_markdown_with_params(results,
                                             markdown_path ="rmarkdown/Summarize_Filtering.Rmd",
                                             dest_path = dest_path,
                                             dest_file_name = dest_file_name,
                                             workdir = workdir
  )
}


