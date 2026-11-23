#' Find a file stored in a package
#' @return The computed result.
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
