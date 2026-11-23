# Add the child-count column to data_ionstar
#
# Problem: the configurations of data_ionstar (config, config_N) name the
# child-count column `nr_children`, but neither $data nor $data_N has it, so
# aggregating the packaged data stopped with "Column `nr_children` not found".
#
# Fix: one row is one peptide of one sample, so every row counts one child.
# Run from the package root after data-raw/fix_deprecated_config.R.

load("data/data_ionstar.rda")
stopifnot(
  data_ionstar$config$nr_children == "nr_children",
  data_ionstar$config_N$nr_children == "nr_children"
)
data_ionstar$data$nr_children <- 1L
data_ionstar$data_N$nr_children <- 1L
save(data_ionstar, file = "data/data_ionstar.rda", compress = "xz")
