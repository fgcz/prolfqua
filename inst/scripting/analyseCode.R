# code snippets to run coverage

library(tidyverse)

# compute test coverage
coverage <- package_coverage(type = "all")

print(coverage, group="functions")
x <- coverage
by = "line"
group = "functions"
df <- tally_coverage(x, by = by)
df <- tidyr::unite(df, "filefunc", filename, functions, remove = FALSE, sep= "~")
group = "filefunc"
percents <- tapply(df$value, df[[group]], FUN = function(x) (sum(x > 0) / length(x)) * 100)

by_coverage <- percents[order(percents,
                              names(percents))]
coveragestats <- data.frame(funcname = names(by_coverage) , coverage = by_coverage)
coveragestats <- coveragestats |> tidyr::separate("funcname",c("file", "functions"), sep="~")
coveragestats |> group_by(functions) |> summarize(n = n()) -> nrfn
nrfn[nrfn$n >= 7,]
coveragestats |> filter(functions == "initialize")

### get function names
functionNames <- data.frame(functionName = ls("package:prolfqua"))

### get documentation
db <- Rd_db("prolfqua")
res <- data.frame(rdfile = names(db),
                  title = unlist(lapply(db, tools:::.Rd_get_metadata, "title")),
                  functionname = unlist(lapply(db, tools:::.Rd_get_metadata, "name")),
                  description = unlist(lapply(db, tools:::.Rd_get_metadata, "description")))
allFun <- data.frame(functionsls = ls("package:prolfqua"))

res  <- full_join(allFun, res, by = c(functionsls= "functionname"))
res <- full_join(coveragestats, res, by = c(funcname = "functionsls"))
View(res)
writexl::write_xlsx(res, path = "inst/scripting/FunctionOverview.xlsx")
