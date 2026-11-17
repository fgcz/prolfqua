# Verification of srr groundwork

**Completed on 2026-11-17.** All four steps and the acceptance checks in [the original plan](../../TODO/prolfqua/Archive/PLAN_srr_groundwork_2026-11-16.md) are complete for prolfqua 1.9.4. Verification used R 4.6.2 on arm64 macOS. Completion of this groundwork does not mean the package is ready for rOpenSci submission; the plan explicitly leaves the design-heavy standards gaps for later work.

## Plan completion

- [x] **Step 1 — deterministic Hurdle fixture:** explicitly remove group A for three suitable proteins and assert those proteins in the missing-group result. Committed as `322f0556`; seeds 1, 2 and 3 were rechecked locally.
- [x] **Step 2 — validation and conditions:** validate the specified entry points, use classed errors, correct duplicated error/warning wording, and test the error and warning paths. Committed as `0c549921`, with follow-up regression fixes for valid GLM ANOVA choices and scalar selectors. All 142 remaining abort/warning sites execute in the final suite.
- [x] **Step 3 — documentation:** document all 107 previously undocumented internal functions, add lifecycle and glossary sections, explain formula/design/contrast coding, record measured runtime scaling, and correct the modelling-method input terminology. Committed as `0adc8a61`, with misleading internal parameter descriptions corrected during review.
- [x] **Step 4 — standards evidence:** reassess all 116 General/Regression standards, add evidence and NA/TODO annotations, update both assessment TSVs and the workspace TODO, and generate the compliance reports. Committed as `1691769f`, with partial claims subsequently changed to TODOs as requested in the review follow-up.

## Final executable checks

| Check | Result |
|---|---|
| `make test` | 1,598 assertions passed; zero failures, warnings or skips |
| Nine moved optional packages absent | Passed; 38 expected skips for unavailable optional functionality |
| Hurdle fixture, seeds 1–3 | 69 assertions passed for each seed; seeds 1–2 had no warnings, seed 3 had four simulated zero-variance warnings |
| `make document` | Passed without warnings |
| `make lint` | Passed |
| Internal documentation audit | No undocumented top-level internal functions remain |
| Package help | Options, Lifecycle and Glossary present in `man/prolfqua-package.Rd` |
| `make vignette V=articles/LowLevelAPI` | Rendered, including the asserted one-factor and interaction calculations |
| `make vignette V=ModellingMethods` | Rendered |
| Standards inventory | All 116 canonical IDs occur exactly once: 33 met, 7 NA and 76 TODO; annotations agree with the assessment TSV |
| `srr::srr_stats_pre_submit()` | All IDs documented; the intended TODOs and below-50% current compliance prevent submission |
| `srr::srr_report()` | Rendered |
| `covr::package_coverage()` | 86.19962% overall; all 142 abort/warning sites executed |
| `make install INSTALL_BUILD_TARGET=` | Installed the source tarball, including vignettes, produced by the final full check |
| Root `make integration SNAKEMAKE='snakemake --forcerun run_integration_test'` | Reran integration tests against the installed package: 119 assertions passed; zero failures, warnings or skips |
| Full `make check` | Zero errors, warnings or notes, including examples, tests and vignette rebuilds |

The committed Hurdle fixture uses seed 7; the complete test suite has no unexpected warnings. The additional seed-3 warnings do not fail its missing-group assertions and are recorded rather than suppressed.

The reduced-library run omitted UpSetR, ComplexHeatmap, circlize, plotly, ggbeeswarm, progress, logistf, lme4 and lmerTest. Packages remained installed in the normal user library; the diagnostic used a separate symlink library and did not move or remove them.

The condition audit is recorded in [condition-coverage.tsv](condition-coverage.tsv). The previously uncovered contrast guard is now exercised by an all-zero contrast. The donor-subset guard was unreachable because its target was selected from existing donor keys; that branch and its unused handler were removed, and covariance-borrowing help now describes the implemented subset selection. Informative-message paths and unique wording are not exhaustively covered, so G5.2/G5.2a/G5.2b remain partial. Those broader message requirements do not change the original plan's completed abort/warning acceptance check.

## Current standards status

Only fully met standards now receive `@srrstats`; partial and unmet requirements retain `@srrstatsTODO`. This follows the user's requested correction of the original plan's partial-tag convention. Current compliance is **20/68 General (29.4%)** and **13/48 Regression (27.1%)**. The plan allows an honest below-50% result; model accessors and retained responses, missing-value policies, convergence controls, broader reference comparisons and recovery tests remain separate follow-up work.

`autotest` mutation testing was not repeated. The recorded investigation exhausted memory twice; targeted input-validation tests and the full package/integration checks above were run instead. The historical traced-test result is not represented as a mutation pass on this revision. See the [assessment](standards-compliance.html) for the recorded limitation.

Review the rendered [assessment](standards-compliance.html), [srr report](srr-report.html), [Under the hood article](../vignettes/articles/LowLevelAPI.html), [Modelling methods vignette](../vignettes/ModellingMethods.html), and downloaded [General](general.html) and [Regression](regression.html) standards.
