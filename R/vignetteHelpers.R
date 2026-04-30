#' find file stored in package
#' @return The computed result.
#' @export
#' @param packagename name of the R package
#' @param file relative path to the file within the package
#' @examples
#' find_package_file("prolfqua","extdata/medata.csv")
#'
find_package_file <- function(packagename, file) {
  src_script <- file.path(find.package(packagename), file)
  if (!file.exists(src_script)) {
    src_script <- file.path(find.package(packagename), "inst", file)
  }
  if (file.exists(src_script)) {
    return(src_script)
  } else {
    return(NULL)
  }
}

#' copy script files and other from a package to workdir
#' @export
#' @keywords internal
#' @examples
#' copied <- script_copy_helper_vec("extdata/metadata.csv", workdir = tempdir())
#' basename(copied)
script_copy_helper_vec <-
  function(runscripts, workdir = getwd(), packagename = "prolfqua") {
    res <- NULL
    for (scripts in runscripts) {
      src_script <- file.path(find.package(packagename), scripts)
      dest_script <- file.path(workdir, basename(scripts))
      message("copy ", src_script, " to ", dest_script)
      if (!file.exists(src_script)) {
        src_script <- file.path(find.package(packagename), "inst", scripts)
        if (!file.exists(src_script)) {
          msg <- sprintf("could not copy script file. %s", dest_script)
          warning(msg, call. = FALSE)
        }
      }
      if (!file.copy(src_script, dest_script, overwrite = TRUE)) {
        msg <- sprintf("could not copy script file. %s to %s", src_script, dest_script)
        warning(msg, call. = FALSE)
      } else {
        res <- c(res, dest_script)
      }
    }
    message(sprintf("your working directory now should contain: %d new files:\n", length(res)))
    return(res)
  }
