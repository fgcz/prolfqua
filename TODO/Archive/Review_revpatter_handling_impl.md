# Review: rev-pattern handling implementation

**Date:** 2026-07-01
**Scope:** Source review of the implementation described in
`prolfqua/TODO/handover_revpattern_handling.md`, covering the relevant `prolfqua`
and `prolfquapp` changes.

No files were changed during the review, and the full test suites were not
rerun. This review is based on reading the implementation, tests, and installed
application entry points.

## Findings

### High: `pattern_decoys = ""` has split semantics between core and app

Core `prolfqua` treats an empty decoy pattern as an explicit opt-in to the
built-in defaults. `AnalysisConfiguration` documents that `""` selects defaults,
`build_contrast_analysis()` drops decoys whenever `pattern_decoys` is non-`NULL`,
and the core test fixture expects decoys to be removed for `pattern_decoys = ""`.

Evidence:

- `prolfqua/R/AnalysisConfiguration.R:57-66`
- `prolfqua/R/build_contrast_analysis.R:105-108`
- `prolfqua/tests/testthat/test-build-contrast-decoy-drop.R:31-44`

`prolfquapp` currently treats the same value as opt-out. `DEAnalyse$build_facade()`
only drops decoys when the pattern is non-empty, and
`test-single-filtering-path.R` explicitly asserts that `pattern_decoys = ""`
keeps decoys in the fit.

Evidence:

- `prolfquapp/R/R6_DEAnalyse.R:169-173`
- `prolfquapp/tests/testthat/test-single-filtering-path.R:78-89`
- `prolfqua/TODO/handover_revpattern_handling.md:123-130`

Impact: the same configuration value can produce different model inputs depending
on whether the caller enters through `prolfqua::build_contrast_analysis()` or
through `prolfquapp::DEAnalyse`. That undermines the stated "one shared detector"
and "single filtering path" contract.

Recommended resolution: choose one public meaning for `pattern_decoys = ""` and
align `AnalysisConfiguration`, `DEAnalyse`, tests, and handover documentation to
that meaning. Given the core detector already safely handles empty patterns as
"defaults only", the least surprising implementation is probably to remove the
`nzchar()` gate in `DEAnalyse$build_facade()` and let non-`NULL` mean opt-in
consistently.

### High/Medium: installed `CMD_DEA_V2.R` still uses the legacy IBAQ pre-filter

The helper implementation was updated so IBAQ works from a copy of the full
`LFQData` rather than pre-filtering through `ProteinAnnotation$clean()`.

Evidence:

- `prolfquapp/R/cmd_helpers.R:461-472`

However, the installed application script still duplicates the older IBAQ code
path and calls:

```r
lfqdataIB <- xd$lfqdata$get_subset(
  xd$protein_annotation$clean(contaminants = GRP2$processing_options$remove_cont)
)
```

Evidence:

- `prolfquapp/inst/application/CMD_DEA_V2.R:219-222`
- `prolfquapp/inst/application/CMD_DEA_V2.R:229-230`

This script is not dead code: the installed `prolfqua_dea.sh` entry point calls
`application/CMD_DEA_V2.R`, and `test-CMD_DEA_V2.R` executes the script but does
not assert decoy/contaminant retention in IBAQ output.

Impact: the main `run_dea()` model/report path is largely on the new single path,
but IBAQ output from the installed DEA app can still follow the retired
`get_subset(clean())` behavior. This contradicts the handover claim that IBAQ was
aligned to the single path.

Recommended resolution: remove the duplicated legacy IBAQ block from
`inst/application/CMD_DEA_V2.R` or delegate output writing to the fixed helper.
Add a script-level regression test that checks decoys/contaminants are not
silently removed before IBAQ export.

### Medium: report text still claims decoys validate the DE false-positive rate

The report text still says decoy sequences are kept because they allow
re-estimating the proportion of falsely identified proteins in the list of
differentially expressed proteins.

Evidence:

- `prolfquapp/vignettes/Grp2Analysis_V2_R6.Rmd:147-151`
- `prolfquapp/tests/testthat/Grp2Analysis_V2_R6.Rmd:150-151`
- `prolfquapp/inst/application/DIANN/_Grp2Analysis_V2.Rmd:136-137`

That statement no longer matches the implementation. Decoys are dropped before
model fitting, have no contrast statistics, and are omitted from volcano and
significant-result sets. The handover explicitly records that a decoy-based
FDR-validation volcano is not available without a separate fit.

Evidence:

- `prolfqua/TODO/handover_revpattern_handling.md:139-140`
- `prolfquapp/R/R6_DEAnalyse.R:169-181`
- `prolfquapp/R/interactive_report_plots.R:94-104`

Impact: users reading the report will infer that decoys remain available for DE
false-positive validation, but the current result tables cannot support that
claim.

Recommended resolution: update the report prose to say decoys are retained in
raw/abundance export for traceability but excluded from the model fit and
therefore absent from DE statistics.

### Medium: contaminant annotation still has an empty-pattern match-all path

The new shared detector code in `prolfqua` protects against empty regex values by
mapping empty, `NULL`, and `"a^"` to safe built-in defaults.

Evidence:

- `prolfqua/R/utilities.R:683-692`
- `prolfqua/R/utilities.R:716-760`

Decoy detection in `ProteinAnnotation` delegates to `prolfqua::is_decoy()`, but
contaminant annotation still uses raw `grepl(self$pattern_contaminants, ...)`.

Evidence:

- `prolfquapp/R/R6_ProteinAnnotation.R:132-144`
- `prolfquapp/R/R6_ProteinAnnotation.R:368-380`

This matters because `preprocess_MSstats_FPDIA()` defaults
`pattern_contaminants = ""` and passes that value into `ProteinAnnotation$new()`.
For a direct `ProteinAnnotation` summary/report path, `grepl("", x)` can label all
proteins as contaminants.

Evidence:

- `prolfquapp/R/preprocess_MSstats.R:89-95`
- `prolfquapp/R/preprocess_MSstats.R:173-182`

Impact: the original empty-pattern class of bug is fixed for decoys but remains
possible for contaminant labels in at least one direct preprocessing path.

Recommended resolution: route contaminant annotation through
`prolfqua::is_contaminant()` as well, or normalize `pattern_contaminants` before
calling `grepl()`. Add a regression test where `pattern_contaminants = ""` does
not mark every protein as contaminant.

## Residual risks / watch items

- `prolfquapp/R/aggregation_IBAQ.R:83-87` still uses an inner join against
  annotation data. It may currently preserve decoys because `ProteinAnnotation`
  seeds `row_annot` from all quant IDs, but future annotation changes could
  reintroduce silent drops. Prefer a test that locks the desired IBAQ row set.
- `ggiraph` is a `Suggests` dependency, but the report chunks call helpers that
  stop when `ggiraph` is unavailable. That may be acceptable if report-rendering
  environments always install Suggests; otherwise the chunks need a graceful
  fallback or conditional execution.
- Direct callers that bypass `prolfquapp::DEAnalyse` and construct facades
  themselves still depend on `config$pattern_decoys` being correctly populated.
  The current app paths stamp it during `DEAnalyse$initialize()`, but this
  remains an integration contract to keep visible.

## Summary

The main direction of the implementation is sound: decoy/contaminant detection
has moved toward shared core helpers, the default app fit path now uses a
targets-only model input, and export joins are mostly preserving raw abundance
rows. The implementation is not yet complete because config semantics are
inconsistent, the installed DEA script still contains a legacy IBAQ filter, and
some reporting/annotation paths still describe or implement the old behavior.

## Resolution (2026-07-01) — prolfquapp `fd20880`

All four findings addressed; prolfquapp suite 261 pass, 0 fail.

- **F1 (High) — resolved.** `DEAnalyse$build_facade()` now gates the targets-only
  drop on `!is.null(pattern_decoys)` (the `nzchar()` check removed), matching
  `prolfqua::build_contrast_analysis`: non-`NULL` (including `""` → defaults) opts
  in; `NULL` is off. The app maps an empty REVpattern to `NULL` upstream, so a
  cleared pattern is off. `test-single-filtering-path.R` split into a `NULL` (off)
  and a `""` (on, defaults) case. Handover updated.
  (`R6_DEAnalyse.R:162-172`, `test-single-filtering-path.R:73-99`)
- **F2 (High/Med) — resolved.** `inst/application/CMD_DEA_V2.R` IBAQ block now uses
  `xd$lfqdata$get_copy()` instead of `get_subset(clean())`, matching the helper.
  (`CMD_DEA_V2.R:219-224`)
- **F3 (Med) — resolved.** Report prose no longer claims decoys validate the DE
  false-positive rate; it states decoys are retained for traceability but excluded
  from the fit (no DE statistics, absent from the volcano). Added a
  decoy-proportion line. (`vignettes/Grp2Analysis_V2_R6.Rmd:147-152`)
- **F4 (Med) — resolved.** `annotate_contaminants()` routes through
  `prolfqua::is_contaminant()`, so an empty `pattern_contaminants` falls back to
  built-in defaults and can never `grepl("", x)` match-all. Regression test added.
  (`R6_ProteinAnnotation.R:366-380`, `test-R6_ProteinAnnotation.R`)

**Residual risks acknowledged (unchanged, tracked in the handover):**
`aggregation_IBAQ.R` still inner-joins the annotation (functional metadata join;
preserves the intended row set today — a lock test is the suggested follow-up);
`ggiraph` remains a Suggests dep (report helpers `stop()` with a clear message if
absent); the `config$pattern_decoys` stamping in `DEAnalyse$initialize` remains the
integration contract for any facade constructed outside `DEAnalyse`.
