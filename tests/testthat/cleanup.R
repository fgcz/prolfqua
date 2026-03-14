if (grepl("/tests/testthat", getwd())) {
  files <- dir()
  if ("cleanup.R" %in% files) {
    toremove <- files[!grepl("*.R$", files)]
    unlink(toremove, recursive = TRUE)
  }
}
