#' Have a hook for the moment to write in different formats
#' @export
#' @examples
#' x<- data.frame(a = 1:3, b=4:6)
#' x
#' assign("lfq_write_format", "xlsx", envir = .GlobalEnv)
#' lfq_write_format
#' #lfq_write_table(x, "test.csv")
#' assign("lfq_write_format", "csv", envir = .GlobalEnv)
#' #lfq_write_table(x, "test.csv")
#'
#' assign("lfq_write_format", "both", envir = .GlobalEnv)
#' #lfq_write_table(x, "test.csv")
lfq_write_table <- function(x, path){
  #lfq_write_format <- match.arg(lfq_write_format)
  if(!exists("lfq_write_format")){
    assign("lfq_write_format", c("xlsx"), envir = .GlobalEnv)
  }
  message("lfq_write_format", lfq_write_format)
  if("csv" %in% lfq_write_format){
    message("writing csv")
    path_csv <- paste0(tools::file_path_sans_ext(path),".csv")
    readr::write_csv(x, path = path_csv)
  }

  if( "xlsx" %in% lfq_write_format ){
    message("writing xlsx")
    path_xlsx <- paste0(tools::file_path_sans_ext(path),".xlsx")
    writexl::write_xlsx(x, path = path_xlsx)
  }
  if( "html" %in% lfq_write_format)
  {
    message("write html")
    path_html <- paste0(tools::file_path_sans_ext(path),".html")
    dt <- DT::datatable(x, filter = 'top', options  = list(pageLength = 50))
    DT::saveWidget(dt, 'foo.html',selfcontained = TRUE)
    file.rename("foo.html", path_html)
  }
}

