#' find file stored in package
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
scriptCopyHelperVec <-
  function(runscripts, workdir = getwd(), packagename = "prolfqua") {
    res <- NULL
    for (scripts in runscripts) {
      src_script <- file.path(find.package(packagename), scripts)
      dest_script <- file.path(workdir, basename(scripts))
      message("copy ", src_script, " to ", dest_script)
      if (!file.exists(src_script)) {
        src_script <- file.path(find.package(packagename), "inst", scripts)
        if (!file.exists(src_script)) {
          warning(paste("could not copy script file.", dest_script, sep = " "))
        }
      }
      if (!file.copy(src_script, dest_script, overwrite = TRUE)) {
        warning(paste("could not copy script file.", src_script, " to ", dest_script, sep = " "))
      } else {
        res <- c(res, dest_script)
      }
    }
    message(paste(
      "your working directory now should contain: ",
      length(res),
      "new files :\n",
      sep = " "
    ))
    return(res)
  }
