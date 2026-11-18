# prolfqua tests

The suite in `tests/testthat` has standard and extended tests. Both run under testthat in the same files.

## Standard tests

- Run with `make test` (or `devtools::test()`); `R CMD check` runs them too
- Use fixed seeds, simulated data and the packaged `data_ionstar` subset
- Optional backends (lme4, lmerTest, logistf, Rfit, limpa, plotting packages) are skipped when not installed

## Extended tests

Extended tests repeat parameter-recovery checks over many seeds, test noise susceptibility for all registered methods, cover nested models on simulated peptide data, grow the number of samples per group, and compare fit times on noiseless and noisy data.

- Switch: environment variable `PROLFQUA_EXTENDED_TESTS=true`; otherwise they are skipped by `skip_if_not_extended()` (`tests/testthat/helper-extended.R`)
- Run locally: `PROLFQUA_EXTENDED_TESTS=true make test`
- CI: `.github/workflows/r.yaml` sets the variable for commits whose message contains `run-extended`
- Runtime: about 5 minutes on a current laptop, most of it in the 50-seed recovery and group-size tests; the standard suite takes about 1.5 minutes
- Memory: under 2 GB
- Platform: any platform that runs the standard suite; timing comparisons use medians of repeated runs
- Downloads: none; data are simulated or packaged
- Artefacts: none; all checks are assertions, nothing needs manual inspection

## Reference outputs

Stored outputs of other software live in `inst/extdata`; the scripts that produce them are in `data-raw`. Do not put data files in `tests/testthat`: `tests/testthat/cleanup.R` deletes every non-`.R` entry there.
