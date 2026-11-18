# Simulated LFQ data together with the true group differences, for parameter-recovery tests.
#
# sim_lfq_data() stores each protein's true group mean in `average_prot_abundance`; the
# sim_lfq_data_*_config() wrappers drop that column, so this helper keeps it. The truth is taken
# from the true means, not from the `fc` argument (whose sign is the reverse of group - Ctrl).

RECOVERY_CONTRASTS <- c(A_vs_Ctrl = "group_A - group_Ctrl", B_vs_Ctrl = "group_B - group_Ctrl")

sim_with_truth <- function(Nprot = 200, N = 4, seed = 1, peptide = FALSE) {
  raw <- .with_seed(seed, sim_lfq_data(Nprot = Nprot, N = N, PEPTIDE = peptide))
  raw$nr_children <- 1
  means <- unique(data.frame(
    protein_Id = paste(raw$proteinID, raw$idtype2, sep = "~"),
    group = raw$group,
    true_mean = raw$average_prot_abundance
  ))
  ctrl <- means[means$group == "Ctrl", c("protein_Id", "true_mean")]
  truth <- merge(means[means$group != "Ctrl", ], ctrl, by = "protein_Id", suffixes = c("", "_ctrl"))
  truth$contrast <- paste0(truth$group, "_vs_Ctrl")
  truth$true_diff <- truth$true_mean - truth$true_mean_ctrl
  setup <- .sim_setup(raw, list(group_ = "group"), peptide = peptide)
  setup$config$is_response_transformed <- TRUE
  list(
    lfqdata = LFQData$new(setup$data, setup$config),
    truth = truth[, c("protein_Id", "contrast", "true_diff")]
  )
}

# Join estimated contrasts with the truth; one row per protein and contrast.
with_truth <- function(contrasts, truth) {
  merge(contrasts[, c("protein_Id", "contrast", "diff")], truth, by = c("protein_Id", "contrast"))
}
