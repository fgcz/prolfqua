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


#' render MQ Summary.
#' @examples
#' run_MQSummary_Rmd(LFQService::sample_analysis,LFQService::skylineconfig, workdir=tempdir(check = FALSE))
#'
run_MQSummary_Rmd <- function(data, config, dest_file_name = "MQSummary.pdf", path =".",  workdir = getwd()){
  markdown_path <-"rmarkdown/MQSummary.Rmd"
  markdown_file <- basename(markdown_path)
  res <- .scriptCopyHelperVec(markdown_path, workdir = workdir)
  rmarkdown::render(res,
                    output_format = rmarkdown::pdf_document(),
                    params=list(data = data, configuration=config$clone(deep=TRUE)), envir = new.env())
  pdf_doc <- paste0(tools::file_path_sans_ext(res),".pdf")

  file.copy(pdf_doc, file.path(path, dest_file_name), overwrite = TRUE)
}
