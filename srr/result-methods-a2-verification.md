# Verification of A2 analysis result methods

Completed implementation and focused standards review on 2026-11-17 for prolfqua 1.9.4, following the [approved A2 design](../../TODO/prolfqua/TODO_srr_block_A2_vcov_fitted_plot.md) and [implementation plan](../../TODO/prolfqua/PLAN_srr_block_A2_vcov_fitted_plot.md). Changes are left uncommitted; no version bump or push is part of this work, following the user's explicit instruction. Verification uses the normal user/system libraries on arm64 macOS with R 4.6.2. Earlier [block A verification](result-methods-verification.html) and [FDR planning verification](fdr-power-verification.html) remain historical records.

## Changes and numerical evidence

- Existing `Model$get_coefficients()` tables now obtain estimates, standard errors and df through the fit accessor. LOD-refitted models use borrowed covariance and df; the existing completeness filter is preserved. Tests compare lm and lmer tables with their summaries, verify borrowed errors directly, and render robust/rank coefficient volcano plots using the normalized columns. The generic table uses the planned two-sided Wald t p-values. `ModelFirth$get_coefficients()` skips failed fits, with a regression test that preserves the remaining subjects. The matrix imputation field now documents bare subject IDs and the actual `estimate_type` flag.
- `vcov()` returns a long table of subject identifiers, `term`, `term2` and `covariance`. Tests compare ordinary-model covariances directly with backend fits, check covariance diagonals against coefficient standard errors across the 20-method registry, and reconstruct limma/vooma covariance from each subject's available cells, weights and design. Failed matrix fits are omitted; combined hurdle results raise the documented component-selection error. Limpa's diagonal errors agreed to numerical precision in the tested fixtures (largest observed difference 5.6e-17); no compensating rescaling was added.
- The default `plot()` displays observed responses and fitted cell means with confidence intervals for up to six subjects selected by the smallest FDR in the first contrast. Explicit subjects and confidence levels are supported. Tests compare lm means and limits with `predict(..., interval = "confidence")` for one- and two-factor designs, mixed-model limits with the Satterthwaite contrast calculation, and limma limits with moderated `topTable()` results. Detection estimates/limits are probabilities; ROPECA facets retain peptide identities. Seven model-family cases are rendered through `ggplot_build()`, and the complete registry exercises default plotting. Hurdle's default remains volcano; an explicit fitted plot directs the user to its components. Volcano and MA alternatives remain available.
- The modelling help and [Modelling methods vignette](../vignettes/ModellingMethods.html) document covariance, plotting controls, internal aggregation and detection encoding, retained responses and output scales. Aggregation and detection encoding have no inverse. Preceding LFQData transformations retain their input columns; no new transformation-history object was introduced.

Two implementation details refine the plan based on executable checks. For multiple factors, `linfct_factors` contains marginal means whose labels do not align with observed factor combinations; the plot uses `linfct_interactions` cell means instead, verified against direct two-factor prediction. For logistf, `fit$terms` contains coefficient names, so its design matrix is constructed from `terms(formula(fit))` and the stored model frame. Both changes preserve the approved fitted-plot behavior.

## Executable verification

| Check | Result |
|---|---|
| Coefficient regression tests | Reproduced the borrowed-error, robust/rank plotting and failed-Firth defects; corrected focused run passes 86 assertions |
| Coefficient/covariance full suite | 2,440 assertions; zero failures, warnings or skips |
| A2 `make test` | 2,517 assertions; zero failures, warnings or skips |
| Latest working-tree `make test` | 2,539 assertions; zero failures, warnings or skips, including the concurrent sample-offset tests |
| Focused facade/result-method file | 725 assertions; zero failures, warnings or skips |
| Reduced-library tests | 2,224 assertions; zero failures/warnings and 46 expected skips |
| `make document` | Generated S3 registrations and help without warnings |
| Full `make check` | Zero errors, warnings or notes, including all examples, tests and rebuilt vignettes |
| Final `make install` | Checked source tarball installed successfully, including vignettes |
| Lint | Checked A2 archive passes `lintr::lint_package()` with the repository’s `.lintr`; see concurrent-work note below |
| Forced root integration test | 119 assertions; zero failures, warnings or skips across all 10 workflow files |
| Fresh `make vignette V=ModellingMethods` | Rendered after deleting the old HTML; fitted and volcano figures visually inspected |
| Standards inventory | All 116 canonical IDs are covered; annotation statuses, both TSVs and Markdown verdicts agree |
| `srr::srr_stats_pre_submit()` | Checked A2 source: 47 met, 7 NA and 62 TODO. Current [raw output](pre-submit.txt), rerun 2026-11-18 after the duplicate G2.4e annotation from the concurrent sample-offset work was removed: 47 met; broader submission gaps remain |
| `srr::srr_report(view = FALSE, roxygenise = FALSE)` | Regenerated and saved [report](srr-report.html) |

The tested A2 source archive has SHA-256 `8e57ad8ee3f5a6f9dbc34b9da60e3d33c10f3a606d4b847725a4c14ddc0d9b5e`. A separate sample-offset feature appeared in `LFQDataTransformer`, its tests and NEWS after that archive was built. Those changes are preserved and are outside this A2 verification. The later working-tree `make lint` reported a complexity finding in the newly expanded transformer class; the checked A2 source passes the same lintr configuration. The full check and reduced-library results above apply to the A2 archive; a subsequent working-tree test run passes all 2,539 assertions, including the sample-offset tests. The concurrent transformer change also repeats G2.4e, already annotated in `tidyMS_data_setup.R`: srr’s raw current-tree count is therefore 48 annotations, while the assessment correctly remains at 47 unique met standards. The A2 archive reports 47. The concurrent source and annotations were not edited as part of A2.

The reduced-library script hides UpSetR, ComplexHeatmap, circlize, plotly, ggbeeswarm, progress, logistf, lme4 and lmerTest through a disposable symlink library, sets `.libPaths()` inside R and confirms those packages cannot be loaded. Skips correspond to unavailable optional functionality. The downstream source scan of prolfquapp and prophosqua found no `plot()` calls on analysis results relying on the previous default; the aggregator's own `$plot()` and explicit contrast plotters are unaffected. No downstream source changes were needed.

## Standards review and remaining scope

The `package-standards` review checked generated S3 registration, internal roxygen, public arguments and returns, runnable examples, optional-dependency behavior, current-version NEWS and vignette semantics. The `stats-standards` review checked covariance/interval sources, subject alignment, multi-factor means, probability transformations, unsupported hurdle operations and numerical backend references. No required correction remains in the implemented scope after the regression and plotting fixes.

RE4.6, RE6.2 and RE4.12 are met. General remains 20/68 (29.4%); Regression rises from 24/48 to 27/48 (56.2%); overall compliance is 47/116 (40.5%), with 7 NA and 62 TODO tags. Covariance-estimation controls, convergence diagnostics and other remaining standards retain their separate TODOs. RE4.12 is assessed against transformations performed by the modelling call, as explicitly settled in the approved design; it does not require a new user-transformation provenance API.

The [compliance assessment](standards-compliance.html) retains the existing submission blockers and recommendations. Whole-package coverage, pkgcheck and remote CI were not rerun for A2. Autotest was not repeated because the recorded attempts exhausted available memory; no successful mutation result is claimed. Historical measurements are not presented as new A2 results.
