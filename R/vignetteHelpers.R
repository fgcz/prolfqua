#' Find a file stored in a package
#'
#' Looks for `file` in the installed package and, for a source package, in its `inst/` folder.
#' @return the path to `file` in the installed package (or in its `inst/` folder of a source package), or
#'   `NULL` if the file does not exist.
#' @export
#' @param packagename name of the R package
#' @param file relative path to the file within the package
#' @examples
#' find_package_file("prolfqua", "CITATION")
#'
find_package_file <- function(packagename, file) {
  src_script <- file.path(find.package(packagename), file)
  if (!file.exists(src_script)) {
    src_script <- file.path(find.package(packagename), "inst", file)
  }
  if (file.exists(src_script)) src_script
}
