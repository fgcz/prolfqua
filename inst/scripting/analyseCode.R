# code snippets to run coverage

library(tidyverse)
library(covr)

# compute test coverage
x <- package_coverage(type = "all", quit=FALSE, clean=FALSE)
x <- coverage
#zero_coverage(x)
by = "line"
group = "functions"
df <- tally_coverage(x, by = by)
df <- tidyr::unite(df, "filefunc", filename, functions, remove = FALSE, sep= "~")
group = "filefunc"
percents <- tapply(df$value, df[[group]], FUN = function(x) (sum(x > 0) / length(x)) * 100)

by_coverage <- percents[order(percents,
                              names(percents))]
coveragestats <- data.frame(funcname = names(by_coverage) , coverage = by_coverage)
coveragestats <- coveragestats |> tidyr::separate("funcname", c("file", "functions"), sep="~")
View(coveragestats)

coveragestats |> group_by(functions) |> summarize(n = n()) -> nrfn



### get documentation
library(prolfqua)

exportedProlfqua <- ls("package:prolfqua")

res <- list()
for(i in exportedProlfqua){
  if(class(get(i)) == "R6ClassGenerator"){
      public_methods <- names(get(i)$public_methods)
      df  <- data.frame(class = i, method = public_methods)
      res[[i]] <- df
  }
}

R6Classes <- dplyr::bind_rows(res)
R6Classes$ItIs <- "R6"
coveragestats <- tibble(coveragestats) |> mutate(file = tools::file_path_sans_ext( basename(file)))
functionsAndTestCoverage <- full_join(R6Classes, coveragestats, by = c(class = "file", method = "functions"), keep = TRUE) |> tibble()

functionsAndTestCoverage <- full_join(tibble(exported = exportedProlfqua), functionsAndTestCoverage, by=c(exported = "functions"), keep=TRUE )

writexl::write_xlsx(functionsAndTestCoverage, path = "inst/scripting/prolfqua_CoverageOverview.xlsx")


db <- tools::Rd_db("prolfqua", dir=".")
docs <- data.frame(rdfile = names(db),
                  title = unlist(lapply(db, tools:::.Rd_get_metadata, "title")),
                  functionname = unlist(lapply(db, tools:::.Rd_get_metadata, "name")),
                  description = unlist(lapply(db, tools:::.Rd_get_metadata, "description")))
allFun <- data.frame(functionsls = ls("package:prolfqua"))
funcdocs  <- full_join(allFun, docs, by = c(functionsls= "functionname"))
writexl::write_xlsx(funcdocs, path = "inst/scripting/prolfqua_DocumentationOverview.xlsx")
