# Simulated LFQ data together with the true group differences, for parameter-recovery tests.
#
# sim_lfq_data() stores each protein's true group mean in `average_prot_abundance` and, for
# peptide data, each peptide's true group mean in `avg_peptide_abd`; the sim_lfq_data_*_config()
# wrappers drop both. The truth is taken from these means, not from the `fc` argument (whose
# sign is the reverse of group - Ctrl).
#
# Truth columns per protein and contrast:
# - true_diff: difference of the true protein means
# - peptide_mean, peptide_median (peptide data only): mean and median over the protein's peptides
#   of their true differences. The peptide simulator draws every peptide's group mean
#   independently, so nested models estimate these rather than true_diff.

RECOVERY_CONTRASTS <- c(A_vs_Ctrl = "group_A - group_Ctrl", B_vs_Ctrl = "group_B - group_Ctrl")

.true_differences <- function(ids, group, value) {
  means <- unique(data.frame(id = ids, group = group, value = value))
  ctrl <- means[means$group == "Ctrl", c("id", "value")]
  diffs <- merge(means[means$group != "Ctrl", ], ctrl, by = "id", suffixes = c("", "_ctrl"))
  diffs$contrast <- paste0(diffs$group, "_vs_Ctrl")
  diffs$diff <- diffs$value - diffs$value_ctrl
  diffs
}

sim_with_truth <- function(Nprot = 200, N = 4, seed = 1, peptide = FALSE, counts = FALSE) {
  raw <- .with_seed(seed, sim_lfq_data(Nprot = Nprot, N = N, PEPTIDE = peptide))
  protein_id <- paste(raw$proteinID, raw$idtype2, sep = "~")
  # counts = TRUE uses the simulated peptide counts as child counts (needed by deqms)
  raw$nr_children <- if (counts) raw$nrPeptides else 1
  proteins <- .true_differences(protein_id, raw$group, raw$average_prot_abundance)
  truth <- data.frame(protein_Id = proteins$id, contrast = proteins$contrast, true_diff = proteins$diff)
  if (peptide) {
    peptides <- .true_differences(paste(protein_id, raw$peptideID, sep = "|"), raw$group, raw$avg_peptide_abd)
    peptides$protein_Id <- sub("\\|.*$", "", peptides$id)
    by_protein <- split(peptides$diff, paste(peptides$protein_Id, peptides$contrast))
    key <- paste(truth$protein_Id, truth$contrast)
    truth$peptide_mean <- vapply(by_protein[key], mean, numeric(1))
    truth$peptide_median <- vapply(by_protein[key], stats::median, numeric(1))
  }
  setup <- .sim_setup(raw, list(group_ = "group"), peptide = peptide)
  setup$config$is_response_transformed <- TRUE
  list(lfqdata = LFQData$new(setup$data, setup$config), truth = truth)
}

# Join estimated contrasts with one truth column; one row per protein and contrast.
with_truth <- function(contrasts, truth, column = "true_diff") {
  truth <- data.frame(protein_Id = truth$protein_Id, contrast = truth$contrast, true_diff = truth[[column]])
  merge(contrasts[, c("protein_Id", "contrast", "diff")], truth, by = c("protein_Id", "contrast"))
}

# Bias, slope of estimate on truth and RMSE, per contrast.
recovery_stats <- function(joined) {
  do.call(
    rbind,
    lapply(split(joined, joined$contrast), function(d) {
      data.frame(
        contrast = d$contrast[1],
        bias = mean(d$diff - d$true_diff),
        slope = unname(stats::coef(stats::lm(diff ~ true_diff, data = d))[2]),
        rmse = sqrt(mean((d$diff - d$true_diff)^2))
      )
    })
  )
}
