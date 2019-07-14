#' Have a hook for the moment to write in different formats
#' @export
#' @examples
#' x<- data.frame(a = 1:3, b=4:6)
#' x
#' assign("lfq_write_format", "xlsx", envir = .GlobalEnv)
#' lfq_write_format
#' lfq_write_table(x, "test.csv")
#' assign("lfq_write_format", "csv", envir = .GlobalEnv
#' lfq_write_table(x, "test.csv")
#'
#' assign("lfq_write_format", "both", envir = .GlobalEnv
#' lfq_write_table(x, "test.csv")
lfq_write_table <- function(x, path){
  #lfq_write_format <- match.arg(lfq_write_format)
  if(!exists("lfq_write_format")){
    assign("lfq_write_format", "xlsx", envir = .GlobalEnv)
  }
  message("lfq_write_format", lfq_write_format)
  if(lfq_write_format == "csv" | lfq_write_format == "both"){
    message("writing csv")
    path <- paste0(tools::file_path_sans_ext(path),".csv")
    readr::write_csv(x, path = path)
  }

  if(lfq_write_format == "xlsx" | lfq_write_format == "both"){
    message("writing xlsx")
    path <- paste0(tools::file_path_sans_ext(path),".xlsx")
    writexl::write_xlsx(x, path = path)
  }
}

