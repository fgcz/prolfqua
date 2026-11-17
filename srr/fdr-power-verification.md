# Verification of FDR sample-size planning

Implementation and focused standards review on 2026-11-17 for prolfqua 1.9.4, following Part A of [the approved plan](/Users/wolski/.claude/plans/yes-please-snuggly-graham.md). Verification uses R 4.6.2 on aarch64 macOS, normal user/system libraries, limma 3.68.5 and ssize.fdr 1.3. The litmusCore design note records the revised ownership decisions; its Part B implementation is separate. Version 1.9.5 and the GitHub push require the approval specified in the plan.

## Implementation

`fdr_sig_level()` computes the Jung/Liu–Hwang per-test planning threshold. Existing per-protein and SD-quantile methods accept it through `sig.level`. `power_fdr_prior_n()` fits an inverse-gamma prior to unweighted, unmoderated pilot variance estimates using `limma::squeezeVar(robust = FALSE)`, then averages ordinary two-sided t-test power over independent variance and normal-effect priors. `LFQDataStats$power_fdr_pilot()` supplies pooled within-group variances and observed degrees of freedom from the configured LFQData groups, independently of the decorator's summary mode.

All prospective sample counts are per group in a balanced two-group design. `delta_se` denotes the standard deviation of true non-null effects, not the standard error of an estimated mean. The user supplies the effect distribution and null-feature proportion. At least three usable features, including two positive variances, are required because the limma prior estimate is otherwise unavailable. Non-finite variance/df estimates and non-positive df are excluded; negative finite variances and malformed controls raise classed errors. Zero variances remain available to limma, with its warnings visible.

Midpoint quadrature uses 256 gamma-precision nodes, increased to 1024 below four prior degrees of freedom, and 512 normal-effect nodes. Infinite prior df reduces to a fixed variance; zero effect SD reduces to a fixed effect. The curve covers every integer n from two through `max_n`. If no point reaches the target, `N` is `NA_integer_` and achieved power is reported at `max_n`.

## Numerical verification

- The documented standardized effect of one, null proportion 0.9, FDR 0.05 and power 0.8 gives alpha 0.00467836257309942 and 29 samples per group.
- A 36-case fixed-effect grid agrees exactly with `ssize.fdr::ssize.twoSamp`: effects 0.5/1/2, null proportions 0.8/0.9/0.99, FDR 0.05/0.1 and power 0.8/0.9.
- A 12-case varying-prior grid agrees within one sample with successful `ssize.twoSampVary` calculations: generating prior df 4/12/40, effect SD 0/0.3 and null proportion 0.9/0.95. The reference receives the same fitted prior parameters. Optional reference tests close their own temporary graphics device and suppress only the reference package's known precision warning; production calculations do not suppress warnings.
- The infinite-df reduction agrees with the entire `stats::power.t.test(strict = TRUE)` curve to tolerance 1e-12. Independent analytic normal mixing at fixed variance exposed inadequate 128-node accuracy; increasing to 512 effect nodes reduced the largest observed absolute power error from 8.84e-4 to 1.39e-4 in that test.
- For prior df two, an independent adaptive integral with analytic normal mixing agrees with the gamma/normal quadrature to absolute error below 2.9e-5 at n = 4, 20 and 60. A prior df of 200 approaches fixed variance without unstable integration.
- A seeded sample of 20,000 pilot variances recovers the generating prior scale within 5% and prior df within 10%. Tests also cover negative-effect symmetry, both rejection tails, unreachable targets, a one-point n grid, filtering and original-input retention, pooled observed df, classed invalid-input errors and the untransformed-data warning/NULL return.

These are numerical checks on the stated model, not empirical guarantees for all pilot distributions. The fixed alpha is calculated at target average power; points below target need not achieve the FDR planning bound. Ordinary t-test degrees of freedom make this sample-size calculation conservative relative to moderation under the same model. Independence or an adequate weak-dependence approximation, equal group variances and the supplied effect/null assumptions remain necessary. The default quadrature is approximate and can affect integer decisions near a power threshold.

## Executable checks

| Check | Result |
|---|---|
| Focused power and warning tests | 183 assertions; zero failures, warnings or skips |
| `make test` | 2,323 assertions; zero failures, warnings or skips |
| `make document` | Generated exports/help without warnings |
| `make lint` | No lints |
| `make check-fast` | Zero errors, warnings or notes, including examples and tests against the built package |
| `make install` | Final source tarball built with all vignettes and installed successfully |
| `make vignette V=QualityControlAndSampleSizeEstimation` | Fresh [HTML vignette](../vignettes/QualityControlAndSampleSizeEstimation.html), with resolved captions and inspected power figure |
| Standards inventory | All 116 canonical IDs occur exactly once; annotations and both assessment TSVs agree |
| `srr::srr_stats_pre_submit()` | 44 met, 7 NA and 65 TODO; existing submission gaps retained in [output](pre-submit.txt) |
| `srr::srr_report(view = FALSE, roxygenise = FALSE)` | Regenerated [report](srr-report.html) |

The rendered pilot example returns 29 samples per group, achieved power 0.8049355, prior variance 0.006989028, prior df 24.91456 and 28 usable features. The render emits optional magick image-cropping notices; figures are still generated and embedded successfully. Full `make check`, cross-package integration, remote CI, whole-package coverage and pkgcheck were not rerun for this focused plan. Autotest was not repeated because the earlier attempts exhausted available memory; no mutation-test pass is claimed.

## Standards review

The `package-standards` review checked the two new exports, generated NAMESPACE and help, runnable examples, R6 method documentation, current-version NEWS, optional test-only dependency and the rendered vignette. The `stats-standards` review checked primary references, scalar/vector/type contracts, statistical assumptions, return semantics and independent numerical tests. No required correction remains in this change after the quadrature and vignette fixes.

The compliance update extends evidence for G1.0, G2.0, G2.1, G5.2 and G5.4 without promoting partial package-wide requirements to met. General remains 20/68 and Regression 24/48; overall 44/116, with 7 NA and 65 TODO. The two missing pkgcheck documentation items, broader CI recommendations and remaining statistical requirements are preserved in the [compliance assessment](standards-compliance.html). They are separate submission work, not completion criteria for this implementation.
