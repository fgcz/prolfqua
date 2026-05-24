# TODO: Contrast-Level ORA and Rank Interfaces

## Goal

Move ORA and rank-list extraction into the `ContrastsInterface` abstraction so downstream report code can write enrichment input files without knowing model-specific score and significance columns.

## Plan

1. Extend `ContrastsInterface` with `get_rank(score = "statistic")` and `get_ora(up = TRUE, FDR_threshold = 0.05, diff_threshold = 1)`.
2. Implement default behavior in `ContrastsTable` for standard contrast schemas using `statistic`, `FDR`, and `diff`.
3. Override the methods in `prolfquasaint::ContrastsSAINTexpress` using `log2_EFCs` and `BFDR`.
4. Update `prolfquapp::DEAReportGenerator` to call the contrast methods and write regular DEA and SAINT ORA/RNK files from their returned tables.
5. Add focused tests in all three packages and reinstall in dependency order with package-local Makefiles.
