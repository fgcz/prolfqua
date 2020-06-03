#' Facade to writing data frame into various file formats.
#'
#' @export
#' @examples
#' x<- data.frame(a = 1:3, b = 4:6 )
#'
#' lfq_write_table(x, "test.csv", lfq_write_format = c( "csv", "xlsx",  "html"))
#'
lfq_write_table <- function(x, path, lfq_write_format = c("xlsx","csv","html"))
{
  lfq_write_format <- match.arg(lfq_write_format, several.ok = TRUE)
  message("lfq_write_format ", paste(lfq_write_format, collapse = ", "))
  if ("csv" %in% lfq_write_format) {
    message("writing csv")
    path_csv <- paste0(tools::file_path_sans_ext(path),".csv")
    readr::write_csv(x, path = path_csv)
  }
  if ( "xlsx" %in% lfq_write_format ) {
    message("writing xlsx")
    path_xlsx <- paste0(tools::file_path_sans_ext(path),".xlsx")
    writexl::write_xlsx(x, path = path_xlsx)
  }
  if ("html" %in% lfq_write_format)
  {
    message("write html")
    path_html <- paste0(tools::file_path_sans_ext(path),".html")
    dt <- DT::datatable(x, filter = 'top', options  = list(pageLength = 50))
    DT::saveWidget(dt, 'foo.html',selfcontained = TRUE)
    file.rename("foo.html", path_html)
  }
}

