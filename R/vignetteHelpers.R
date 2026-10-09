#' find file stored in package
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

#' copy script files and other from a package to workdir
#' @param runscripts paths of the files relative to the package root (or \code{inst/})
#' @param workdir directory to copy the files to
#' @param packagename package to copy the files from
#' @param overwrite overwrite existing files in \code{workdir}, default TRUE
#' @return paths of the copied files
#' @export
#' @keywords internal
#' @examples
#' copied <- script_copy_helper_vec("CITATION", workdir = tempdir())
#' basename(copied)
script_copy_helper_vec <-
  function(runscripts, workdir = getwd(), packagename = "prolfqua", overwrite = TRUE) {
    res <- NULL
    for (scripts in runscripts) {
      src_script <- find_package_file(packagename, scripts)
      dest_script <- file.path(workdir, basename(scripts))
      message("copy ", src_script, " to ", dest_script)
      if (is.null(src_script)) {
        warning(sprintf("could not copy script file. %s", dest_script), call. = FALSE)
      } else if (!file.copy(src_script, dest_script, overwrite = overwrite)) {
        warning(sprintf("could not copy script file. %s to %s", src_script, dest_script), call. = FALSE)
      } else {
        res <- c(res, dest_script)
      }
    }
    message(sprintf("your working directory now should contain: %d new files:\n", length(res)))
    return(res)
  }
