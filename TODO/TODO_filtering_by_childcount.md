# TODO: generalize child-count filtering in `LFQData`

> `filter_proteins_by_peptide_count()` does not fit `LFQData`'s generic,
> any-hierarchy-level nature, and its threshold lives in the wrong place
> (`config$min_peptides_protein`). Reconsider both.

## Context

`prolfquapp` is adding a minimum-peptides-per-parent-protein filter and
deliberately does **not** reuse this prolfqua function (see
`prolfquapp/TODO/TODO_expose_nr_peptides.md`). While investigating, two problems
with the existing prolfqua API surfaced. This TODO records them for a separate
prolfqua-side cleanup; it is not blocking the prolfquapp work.

## Requirements

- The child-count filter should match `LFQData`'s design: it works at **any
  hierarchy level** (protein→peptide, protein→precursor, peptide→fragment, …),
  parameterized by an explicit **parent key** and **child key**, not tied to the
  modelling `hierarchy_depth`.
- The threshold should be an **explicit argument**, not read implicitly from a
  config field.
- Backwards compatibility: the current no-arg `LFQData$filter_proteins_by_peptide_count()`
  is used by the `prolfquabenchmark` vignettes and the prolfqua test/vignette —
  their behavior must be preserved (or migrated deliberately).

## Current state (verified in code)

- Exported `filter_proteins_by_peptide_count(pdata, config)`
  (`R/LFQData.R:549-561`) delegates to `nr_B_in_A(pdata, config)`
  (`R/LFQData.R:520-529`): `level_a = config$hierarchy_keys_depth()`,
  `level_b = hierarchy_keys()[length(level_a) + 1]`, then keeps rows with
  `count >= config$min_peptides_protein` (`R/LFQData.R:554`).
- **Depth coupling:** it counts the hierarchy level *just below the modelling
  depth*. For a peptide-level (`_PEPTIDE`) config with `hierarchy_depth = 2`
  there is no level below → `nr_B_in_A` warns "here is no B in A" and returns
  `NULL` → the filter is a **no-op**. Confirmed empirically in the review.
- **Threshold source:** `config$min_peptides_protein` defaults to `2`
  (`R/AnalysisConfiguration.R:79-80`) — an implicit config field rather than an
  argument, so callers can't ask for a different level/threshold without mutating
  config.
- **`LFQData$filter_proteins_by_peptide_count()`** (`R/LFQData.R:181-185`) is the
  no-arg method wrapper; used by `prolfquabenchmark` vignettes (multiple),
  `prolfqua` `test-LFQData.R:21`, and `vignettes/SimulateData.Rmd`.

## Design (sketch)

- Rework the core to accept **explicit `parent_key`, `child_key`, and `threshold`**
  arguments, independent of `hierarchy_depth`; count distinct `child_key` per
  `parent_key`, keep parents at/above `threshold`.
- Keep a thin backwards-compatible `LFQData$filter_proteins_by_peptide_count()`
  that defaults `parent_key = hierarchy_keys()[1]`, `child_key = hierarchy_keys()[2]`,
  `threshold = config$min_peptides_protein`, so existing vignette/test call sites
  keep working.
- Decide the fate of `config$min_peptides_protein`: keep it as the default
  source for the no-arg method, or deprecate in favor of an explicit argument.

## Implementation plan (later)

- [ ] Add explicit-key/threshold arguments to `filter_proteins_by_peptide_count()`
      (and/or `nr_B_in_A`), independent of `hierarchy_depth`.
- [ ] Make the no-arg `LFQData$` method delegate with backwards-compatible
      defaults; verify `prolfquabenchmark` vignettes and `test-LFQData.R` still pass.
- [ ] Decide whether to keep or deprecate `config$min_peptides_protein`.
- [ ] Tests: filtering works at protein→peptide *and* peptide-level/nested shapes
      (the `_PEPTIDE` case that currently no-ops), at explicit thresholds.

## Open questions

- Should the explicit-key filter be the primitive and the config field become
  just a default, or should `min_peptides_protein` be removed entirely?
- Does `prolfquapp`'s reader-local helper eventually get promoted into this
  prolfqua primitive once it is generic, so there is one implementation?
