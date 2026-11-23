# Reproduce the runtime table from the prolfqua repository root.
# Run: Rscript srr/benchmark-runtime.R
# Mixed-model convergence warnings are retained for inspection.
devtools::load_all(quiet = TRUE)
options(prolfqua.verbose = FALSE, prolfqua.progress = function(...) invisible(NULL))
rows <- list()
for (n in c(500L, 1000L, 2000L)) {
  for (method in c("lm", "limma", "lmer_nested")) {
    set.seed(42)
    nested <- method == "lmer_nested"
    data <- expand.grid(
      protein_Id = paste0("p", seq_len(n)),
      sampleName = paste0("s", 1:12),
      peptide_Id = paste0("pep", seq_len(if (nested) 3 else 1)),
      stringsAsFactors = FALSE
    )
    data$sample <- data$sampleName
    data$group_ <- rep(c("A", "B", "Ctrl"), each = 4)[match(data$sample, paste0("s", 1:12))]
    data$abundance <- 20 +
      rnorm(n)[match(data$protein_Id, paste0("p", seq_len(n)))] +
      c(A = 1, B = -1, Ctrl = 0)[data$group_] +
      rnorm(n * 3, sd = 0.5)[match(
        paste(data$protein_Id, data$peptide_Id),
        paste(rep(paste0("p", seq_len(n)), 3), rep(paste0("pep", 1:3), each = n))
      )] +
      rnorm(12, sd = 0.2)[match(data$sample, paste0("s", 1:12))] +
      rnorm(nrow(data), sd = 0.2)
    data$qValue <- 0
    data$nr_children <- if (nested) 1 else 3
    cfg <- AnalysisConfiguration$new()
    cfg$hierarchy <- if (nested) {
      list(protein_Id = "protein_Id", peptide_Id = "peptide_Id")
    } else {
      list(protein_Id = "protein_Id")
    }
    cfg$factors <- list(group_ = "group_")
    cfg$file_name <- "sample"
    cfg$is_response_transformed <- TRUE
    cfg$set_response("abundance")
    lfq <- LFQData$new(data, cfg)
    gc()
    timing <- system.time(
      fit <- build_contrast_analysis(lfq, "~ group_", c(A_vs_Ctrl = "group_A - group_Ctrl"), method = method)
    )
    rows[[length(rows) + 1L]] <- data.frame(
      method = method,
      proteins = n,
      observations = nrow(data),
      elapsed = unname(timing[["elapsed"]])
    )
    write.table(do.call(rbind, rows), "srr/runtime.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
writeLines(c(R.version.string, capture.output(Sys.info()), capture.output(sessionInfo())), "srr/runtime-session.txt")
