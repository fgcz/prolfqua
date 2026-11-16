# Verification of srr groundwork

Verified locally on 2026-11-16 with R 4.6.2 on arm64 macOS.

| Check | Result |
|---|---|
| `make test` | 1,566 assertions passed; zero failures, warnings or skips |
| Nine moved optional packages absent | Passed; 38 expected skips |
| Hurdle fixture, seeds 1–3 | Passed; expected simulated zero-variance warnings observed |
| `make document` | Passed without warnings |
| `make lint` | Passed; existing R6-class exclusions followed their shifted line numbers |
| Internal documentation audit | All 107 missing top-level blocks added; no undocumented top-level internal functions remain |
| Parsed-code comparison | Documentation and standards annotations leave executable R syntax unchanged |
| `make vignette V=articles/LowLevelAPI` | Rendered, including the asserted one-factor and interaction calculations |
| `make vignette V=ModellingMethods` | Rendered |
| Standards inventory | All 116 canonical IDs occur exactly once: 88 evidence annotations, 8 NA, 20 TODO |
| `srr::srr_stats_pre_submit()` | All IDs documented; 20 intended TODOs prevent submission |
| `srr::srr_report()` | Rendered |
| `covr::package_coverage()` | 86.08929% overall; 140 of 142 abort/warning sites executed |
| Root `make integration` | Original target was cached; explicitly forced `run_integration_test` and passed 119 assertions with zero failures, warnings or skips |
| Full `make check` | Passed: zero errors, warnings or notes, including examples, tests and vignette rebuilds |

The reduced-library run omitted UpSetR, ComplexHeatmap, circlize, plotly, ggbeeswarm, progress, logistf, lme4 and lmerTest. Packages remained installed in the normal user library; the diagnostic used a separate symlink library and did not move or remove them.

The full check exposed a missing R6 dependency declaration and namespace import. Both are now explicit: R6 is in Imports and roxygen generates the R6Class namespace import. The final full check passes with no errors, warnings or notes.

Two internal invariant guards remain unexercised: the no-donor-subset guard in `compute_borrowed_variance()` after choosing a subset from the existing donor names, and the coefficient-order guard in `compute_contrast()` after indexing each structure by the same coefficient names. Their absence from execution coverage is explicit in [condition-coverage.tsv](condition-coverage.tsv); no artificial test forces inconsistent base indexing just to increase a counter. Informative-message paths and unique wording are not exhaustively covered, so G5.2/G5.2a/G5.2b remain partial.

Generated HTML was checked for missing local image/script resources. Review the [assessment](standards-compliance.html), [srr report](srr-report.html), [Under the hood article](../vignettes/articles/LowLevelAPI.html), [Modelling methods vignette](../vignettes/ModellingMethods.html), and downloaded [General](general.html) and [Regression](regression.html) standards.
