# Run from the package root with the normal user/system libraries.
# A disposable symlink library hides the nine moved optional dependencies.
excluded <- c('UpSetR', 'ComplexHeatmap', 'circlize', 'plotly', 'ggbeeswarm', 'progress', 'logistf', 'lme4', 'lmerTest')
minlib <- tempfile('prolfqua-minlib-')
dir.create(minlib)
for (library in .libPaths()) {
  packages <- list.dirs(library, recursive = FALSE, full.names = TRUE)
  for (path in packages) {
    name <- basename(path)
    target <- file.path(minlib, name)
    if (!name %in% excluded && !file.exists(target)) stopifnot(file.symlink(path, target))
  }
}
.libPaths(c(minlib, .Library))
stopifnot(!any(vapply(excluded, requireNamespace, logical(1), quietly = TRUE)))
cat('Excluded packages verified absent:', paste(excluded, collapse = ', '), '\n')
results <- devtools::test('.', reporter = 'summary', stop_on_failure = TRUE)
counts <- as.data.frame(results)
print(colSums(counts[c('failed', 'error', 'warning', 'skipped', 'passed')]))
stopifnot(sum(counts$failed) == 0, !any(counts$error), sum(counts$warning) == 0)
unlink(minlib, recursive = TRUE)
