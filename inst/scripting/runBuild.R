#!/bin/bash
#Sys.setenv(R_LIBS_SITE="/scratch/PROLFQUA/r-site-library/")
if(dir.exists("./test_build/")){
  unlink("./test_build",recursive = TRUE, force = TRUE )
}

dir.create("./test_build/")
setwd("./test_build/")

message(">>> cloning repository")
retval = system2("git", args = c("clone", "https://github.com/wolski/prolfqua"))
if(retval != 0){
  stop("Could not clone!")
}

message(">>> building package")
retval = system2("R", args = c("CMD","build", "--log", "prolfqua"))
if(retval != 0){
  stop("Could not build prolfqua!")
}

message(">>> running package check")
prolfqtar <- dir(".",pattern = "prolfqua_0.*.tar.gz")
retval = system2("R", args = c("CMD","check", prolfqtar))
if(retval != 0){
  stop("package check failed")
}

message(">>> running biocheck")
BiocCheck::BiocCheck(prolfqtar)
message(">>> running packagedown")

pkgdown::build_site(pkg = 'prolfqua')

setwd("..")

