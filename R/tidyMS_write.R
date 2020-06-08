#' Facade to writing data frame into various file formats.
#'
#' @export
#' @keywords internal
#' @param x data.frame
#' @param path
#' @param name file name
#' @param lfq_write_format formats to write to.
#' @return list with paths to files written
#' @examples
#' x<- data.frame(a = 1:3, b = 4:6 )
#'
#' # lfq_write_table(x, "." , "test", lfq_write_format = c( "csv", "xlsx",  "html"))
#'
lfq_write_table <- function(x, path, name, lfq_write_format = c("xlsx","csv","html"))
{
  lfq_write_format <- match.arg(lfq_write_format, several.ok = TRUE)
  file_paths <- list()
  if ("csv" %in% lfq_write_format) {
    message("writing csv")
    path_csv <- file.path(path,paste0(name,".csv"))
    file_paths$csv <- path_csv
    readr::write_csv(x, path = path_csv)
  }
  if ( "xlsx" %in% lfq_write_format ) {
    message("writing xlsx")
    path_xlsx <- file.path(path,paste0(name,".xlsx"))
    file_paths$xlsx <- path_xlsx

    writexl::write_xlsx(x, path = path_xlsx)
  }
  if ("html" %in% lfq_write_format)
  {
    message("write html")
    path_html <- file.path(path,paste0(name,".html"))
    dt <- DT::datatable(x, filter = 'top', options  = list(pageLength = 50))
    DT::saveWidget(dt, 'foo.html',selfcontained = TRUE)
    file.rename("foo.html", path_html)
    file_paths$html <- path_html
  }
  invisible(file_paths)
}

