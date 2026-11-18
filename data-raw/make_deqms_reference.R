# Stored reference output of DEqMS::spectraCounteBayes() (srr G5.4b).
#
# Used by tests/testthat/test-reference-implementations.R, which compares prolfqua's
# `deqms` method against it. DEqMS is not a prolfqua dependency, so the reference is
# computed once here and stored in inst/extdata/deqms_reference.rds.
#
# The data follow DEqMS's premise: the residual variance falls with the peptide count
# and scatters around that trend (scaled inverse chi-square with 6 df). sim_lfq_data()
# is not used because it gives every protein the same residual SD.
#
# Run from the package root: Rscript data-raw/make_deqms_reference.R

pkgload::load_all(quiet = TRUE)

simulate_count_variance <- function(n_prot = 300, n_rep = 4, prior_df = 6) {
  groups <- c("Ctrl", "A", "B")
  count <- stats::rgeom(n_prot, 0.3) + 1
  prior_sd <- 0.8 * count^(-0.5)
  variance <- prior_sd^2 * prior_df / stats::rchisq(n_prot, prior_df)
  effect <- matrix(0, n_prot, length(groups), dimnames = list(NULL, groups))
  effect[1:30, "A"] <- 1
  effect[31:60, "B"] <- -1
  base <- stats::rnorm(n_prot, 20, 2)
  design <- expand.grid(replicate = seq_len(n_rep), group = groups, protein = seq_len(n_prot))
  data.frame(
    proteinID = sprintf("P%03d", design$protein),
    idtype2 = "sim",
    nrPeptides = count[design$protein],
    sample = paste0(design$group, "_V", design$replicate),
    group = design$group,
    abundance = base[design$protein] +
      effect[cbind(design$protein, match(design$group, groups))] +
      stats::rnorm(nrow(design), 0, sqrt(variance[design$protein]))
  )
}

raw <- .with_seed(2026, simulate_count_variance())
wide <- tidyr::pivot_wider(
  raw,
  id_cols = c("proteinID", "idtype2", "nrPeptides"),
  names_from = "sample",
  values_from = "abundance"
)
samples <- setdiff(names(wide), c("proteinID", "idtype2", "nrPeptides"))
expression <- as.matrix(wide[, samples])
rownames(expression) <- paste(wide$proteinID, wide$idtype2, sep = "~")
group <- factor(sub("_V[0-9]+$", "", samples))
design <- stats::model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- limma::lmFit(expression, design)
contrast_matrix <- limma::makeContrasts(A_vs_Ctrl = A - Ctrl, B_vs_Ctrl = B - Ctrl, levels = design)
fit <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))
fit$count <- wide$nrPeptides
fit <- DEqMS::spectraCounteBayes(fit, coef_col = 1:2)

reference <- do.call(
  rbind,
  lapply(colnames(contrast_matrix), function(contrast) {
    data.frame(
      protein_Id = rownames(expression),
      contrast = contrast,
      logFC = unname(fit$coefficients[, contrast]),
      sca.t = unname(fit$sca.t[, contrast]),
      sca.p = unname(fit$sca.p[, contrast]),
      sca.df = unname(fit$sca.dfprior + fit$df.residual)
    )
  })
)

saveRDS(
  list(
    raw = raw,
    reference = reference,
    d0 = fit$sca.dfprior,
    versions = c(
      DEqMS = as.character(utils::packageVersion("DEqMS")),
      limma = as.character(utils::packageVersion("limma"))
    )
  ),
  file = "inst/extdata/deqms_reference.rds"
)
