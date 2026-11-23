# Compare prolfqua's ropeca_nested with PECA, the original implementation (Suomi and Elo 2017)
#
# Question: does `ropeca_nested` reproduce PECA's peptide-median p-values, and where does it differ?
# Data: data_ionstar (filtered), dilutions a and e, peptides observed in all 8 samples, log2 intensities.
# PECA: PECA_df(test = "t", type = "median", log = FALSE): ordinary two-sample peptide t-tests, median t per
# protein, p-value of the median t, beta correction pbeta(p, n/2 + 0.5, n/2 + 0.5) with n = number of peptides.
# prolfqua: ropeca_nested on the same peptides: per-peptide lm contrasts, median of the signed peptide p-values,
# the same beta correction with n capped at max.n = 10.
#
# Run from the package root: Rscript data-raw/compare_peca.R
# Needs the Bioconductor package PECA. Results are quoted on ?ContrastsROPECANestedFacade.

devtools::load_all(quiet = TRUE)
options(prolfqua.verbose = FALSE)
cat("prolfqua", as.character(packageVersion("prolfqua")), "PECA", as.character(packageVersion("PECA")), "\n")

ionstar <- prolfqua_data("data_ionstar")$filtered()
data <- dplyr::filter(ionstar$data, .data$dilution. %in% c("a", "e"))
lfq <- LFQData$new(data, ionstar$config)
lfq <- lfq$get_Transformer()$log2()$lfq

# complete peptides only: PECA's row t-tests need every value
wide <- lfq$data_wide(as.matrix = TRUE)
complete <- stats::complete.cases(wide$data)
lfq$set_data(dplyr::semi_join(
  lfq$data_long(),
  wide$rowdata[complete, ],
  by = lfq$hierarchy_keys()
))
wide <- lfq$data_wide(as.matrix = TRUE)
cat("peptides", nrow(wide$data), "proteins", length(unique(wide$rowdata$protein_Id)), "\n")

# prolfqua
ropeca <- build_contrast_analysis(lfq, "~ dilution.", c(e_vs_a = "dilution.e - dilution.a"), "ropeca_nested")
ropeca <- ropeca$get_contrasts()

# PECA, group 1 = e, group 2 = a, so the log ratio is e - a
samples <- wide$annotation
peca_input <- data.frame(protein_Id = wide$rowdata$protein_Id, wide$data, check.names = FALSE)
peca <- PECA::PECA_df(
  peca_input,
  id = "protein_Id",
  samplenames1 = samples$sampleName[samples$dilution. == "e"],
  samplenames2 = samples$sampleName[samples$dilution. == "a"],
  log = FALSE,
  test = "t",
  type = "median"
)
peca$protein_Id <- rownames(peca)

both <- dplyr::inner_join(ropeca, peca, by = "protein_Id", relationship = "one-to-one")
both$group <- dplyr::case_when(
  both$n > 10 ~ "n > 10",
  both$n %% 2 == 1 ~ "odd n <= 10",
  .default = "even n <= 10"
)
summary_table <- both |>
  dplyr::summarise(
    proteins = dplyr::n(),
    max_abs_diff_fc = max(abs(.data$diff - .data$slr)),
    max_abs_diff_p = max(abs(.data$p.value - .data$p)),
    median_ratio_p = stats::median(.data$p.value / .data$p),
    spearman_p = stats::cor(.data$p.value, .data$p, method = "spearman"),
    .by = "group"
  ) |>
  dplyr::arrange(.data$group)
print(as.data.frame(summary_table), digits = 3)
cat(
  "all proteins: spearman p",
  signif(stats::cor(both$p.value, both$p, method = "spearman"), 4),
  "; max |fold change difference|",
  signif(max(abs(both$diff - both$slr)), 3),
  "\n"
)
