# TODO status audit, 2026-07-06

Scope: active top-level `TODO/` files plus inline source TODO markers. Files under
`TODO/Archive/` were treated as historical unless an active file still points to
them. No source code was changed.

Update after cleanup: the addressed files identified below were moved to
`TODO/Archive/`, and `TODO/DONE_summary.md` records the archive decision.

Verification used:

- `rg -n "TODO|todo|FIXME|XXX"` over active source/test/vignette/man/TODO files.
- Focused static checks of the referenced implementations.
- `make check-bioc`.
- Focused tests for rfit-impute, vectorized contrasts, aggregation, and limma.
- Manual reproduction of the vectorized rank-deficient contrast bug.

## Already addressed or mostly resolved

### `TODO/Archive/TODO_rfit_impute.md`

Status: addressed and archived.

Evidence:

- `ContrastsRfitImputeFacade` exists and calls `strategy_rfit()` plus
  `build_model_impute(..., borrow_method = "vcov", on_misalign = "fail")`.
- The facade registry contains `rfit_impute`.
- `tests/testthat/test-ContrastsRfitImpute.R` covers the wrapper, builder,
  facade schema, `get_missing()`, LOD-created ties, null-model coefficient names,
  dispatch, and registry.
- Focused test result: 49 pass, 0 fail, 32 warnings from expected difficult
  missing/rfit fixtures.

Archive note: recorded in `TODO/DONE_summary.md`.

### `TODO/Archive/TODO_volcano_plotly_or_ggplot.md`

Status: addressed by decision and archived; no implementation work is pending in
this note.

The file explicitly records the decision to leave both volcano implementations
as-is. It is a decision record, not an active TODO.

Archive note: recorded in `TODO/DONE_summary.md`.

### `TODO/Archive/Review_revpatter_handling_impl.md`

Status: addressed and archived.

The review findings are followed by a resolution section saying all four findings
were fixed in prolfquapp `fd20880`. Current source also shows the empty-pattern
decoy semantics and contaminant detector paths aligned with that resolution.

Archive note: recorded in `TODO/DONE_summary.md`.

### Revpattern handover/design files

Files:

- `TODO/Archive/TODO_revpattern_handling.md`
- `TODO/Archive/handover_revpattern_handling.md`

Status: mostly addressed and archived as historical context. Remaining
watch-items are future choices, not the original root issue.

Evidence:

- prolfqua core has `is_decoy()`, `is_contaminant()`, `effective_*_pattern()`,
  `pattern_decoys`, `pattern_contaminants`, `remove_decoys()`, QC proportions,
  and the target-only fit gate.
- prolfquapp has pattern stamping, target-only fit, single-filtering-path tests,
  contaminant annotation tests, and report text saying contaminants are flagged
  and drawn as triangles.
- The older design note still says contaminant figure labeling was postponed,
  but the handover and current prolfquapp code show ggiraph report plots were
  later added.

Remaining items are deferred/future choices rather than the original root issue:
nested decoy-proportion, possible future contaminant removal mode, separate
decoy FDR-validation fit, and DIA-NN quant-side contaminant detection if needed.

Archive note: recorded in `TODO/DONE_summary.md`; create separate precise TODOs
only for deferred items that are still wanted.

## Still open or stale

### `TODO/TODO_vectorized_contrast_na_coef.md` (untracked)

Status: open; confirmed still reproducible.

Evidence:

- `R/tidyMS_contrasts.R` still uses a signed sum over NA-coefficient contrast
  weights:
  `linfct[, na_coefs, drop = FALSE] %*% rep(1, sum(na_coefs)) != 0`.
- Manual reproduction gives loop path `estimate = NA, p.value = NA` and
  vectorized path `estimate = 0, p.value = NaN`.
- The existing vectorized test file passes only because the simulated
  rank-deficient case is skipped: "No models with NA coefficients in simulated
  data".

Suggested next step: add the hand-built aliased `lm()` regression test, then fix
the invalid-row detection to check any nonzero weight on an NA coefficient.

### `TODO/Archive/TODO_BiocCheck.md`

Status: archived by maintainer decision.

Current `make check-bioc` result:

```text
1 ERROR | 1 WARNING | 8 NOTES
```

Notable current items:

- ERROR: package still needs to be added to Watched Tags in the Bioconductor
  Support Site profile.
- WARNING: `.Deprecated()` / `.Defunct()` usage found twice.
- NOTES: `cat`/`print`, redundant `stop`/`warn*` in signal conditions,
  `suppressWarnings()`, 39 functions over 50 lines, missing runnable examples
  for `ContrastsFacadeBase.Rd`, `list_facades.Rd`, `unregister_facade.Rd`,
  Bioc-Devel subscription cannot be determined, and formatting notes remain.

Archive note: the support-site watched-tag item will be handled outside this
repo. Other BiocCheck cleanup can be reopened as separate code/doc TODOs if
needed.

### `TODO/todo_long_functions.ms`

Status: open and stale.

BiocCheck currently reports 39 functions over 50 lines. The file says 37 and
also warns that it should be regenerated before the next cleanup pass. A
lightweight parse check also found 39 functions over 50 lines.

Suggested cleanup: regenerate the inventory from the current tree before using
the per-function priorities.

### `TODO/TODO_filtering_by_childcount.md`

Status: open.

Evidence:

- `LFQData$filter_proteins_by_peptide_count()` still takes no arguments and
  delegates to `filter_proteins_by_peptide_count(pdata, config)`.
- `nr_B_in_A()` still derives `level_a` from `config$hierarchy_keys_depth()` and
  `level_b` from the next hierarchy level.
- The threshold still comes from `config$min_peptides_protein`.

The explicit `parent_key`, `child_key`, and `threshold` API has not landed.

### `TODO/TODO_review_nr_of_children.md`

Status: partially addressed, still open.

Addressed:

- Regression tests now pin the intended `nr_obs_sample()` and
  `nr_children_experiment()` semantics.
- `nr_children_experiment()` uses an explicit intermediate
  `nr_children_per_sample`, making the old "maxes the raw column" misread less
  likely.

Still open:

- `nr_features_experiment()` still defaults `name_nr_child = "nr_child_exp"`,
  colliding with `nr_children_experiment()`.
- There is still no consolidated documented taxonomy or shared family wrapper
  for the child-count functions.

### `TODO/TODO.md`

Status by section:

- Limma `block` / `duplicateCorrelation`: open. No `block` or `correlation`
  support is present in `strategy_limma()` or `lmFit()` calls.
- Limma fitting `method = "robust"`: open. `strategy_limma()` has `robust` for
  `eBayes()`, but no `lmFit(method = "robust")` path.
- Test coverage expansion: partly addressed but not closable as written. The
  named test files exist and many newer backend tests exist; the remaining
  request is broad and should be turned into concrete edge-case bullets.
- Code review minor items: partly addressed. `Contrasts_from_factors.Rmd` is no
  longer a stub, and `intensity_summary_by_hkeys()` no longer exists in the tree.
  The generic "missing `@param` documentation" item still needs a current
  documentation audit.
- R6 encapsulation follow-up: open. `AnalysisConfiguration` has no
  `add_hierarchy()` / `add_factor()` methods, and `Model` / `Contrasts` still
  expose public `model_df`, `model_name`, and `subject_id` fields.

### Inline source TODOs / placeholders

Open:

- `R/LFQData.R`: roxygen has `@param is_pep todo`, but `initialize()` has no
  `is_pep` argument. The stale parameter is still generated into `man/LFQData.Rd`.
- `R/tidyMS_data_setup.R`: `# TODO add better warning....` is still present.
- `R/Contrasts.R`: `# TODO (goes into calling code)` is still present and needs
  triage or removal.

Stale reference:

- `R/ContrastsLimma.R` roxygen points to `TODO/TODO_limma_voom_integration.md`,
  but that file now lives in `TODO/Archive/`.

## Suggested cleanup order

1. Refresh `TODO/todo_long_functions.ms` from the current BiocCheck output.
2. Split `TODO/TODO.md` into concrete remaining items, moving the already
   addressed vignette/function-removal notes into `DONE_summary.md`.
3. Fix the untracked vectorized contrast TODO as a normal bug: failing test
   first, root-cause fix in `R/tidyMS_contrasts.R`, then `NEWS.md`.
