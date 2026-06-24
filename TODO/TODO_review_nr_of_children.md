# TODO: Review the "number of children" counting family

## Goal

Consolidate and clarify the family of functions that count "how many B per A"
(peptides per protein, precursors per peptide, observed children per protein
per sample, …). The intuition prompting this review: there are **too many
flavours** of essentially the same operation, with overlapping names, an
outright name collision, and silently different semantics. A caller can easily
pick the wrong one and feed a subtly wrong count into fitting weights or DEqMS
variance moderation.

This is a **review / investigation** task, not an approved refactor. Decide the
target taxonomy first (see "Open questions"), then plan the consolidation
separately.

## Why now

While triaging the code-review item *"`nr_children_experiment` maxes the
un-aggregated column"* (a false positive — see
`../TODO_prolfqua_prolfquapp_code_review.md`), the only reason it *looked* like
a bug was the confusing structure of this function family: a per-sample **sum**
written back into the source column name, a sibling that **distinct-counts**
under a near-identical name, and three more functions computing
"peptides per protein" three different ways in three different files. Correct
code that an experienced reader misreads as broken is itself a defect.

## Inventory (the current flavours)

### `R/tidyMS_aggregation.R`
- `nr_obs_sample(data, response, hierarchy_keys_depth, file_name, nr_children_col, new_child = nr_children_col)`
  — groups observed rows (response not NA) by hierarchy×sample and **sums an
  existing `nr_children` count column**. Despite "obs" in the name it does not
  count rows; it sums a pre-existing count. Output column defaults to the
  *source* name (`new_child = nr_children_col`).
- `nr_children_experiment(..., name_nr_child = "nr_child_exp")` — `max()` over
  samples of `nr_obs_sample`. **Observation-aware.** Feeds DEqMS.
- `nr_features_experiment(data, hierarchy_keys, hierarchy_keys_depth, name_nr_child = "nr_child_exp")`
  — **distinct count** of child features per hierarchy unit. **NOT
  observation-aware** (ignores response). ⚠️ **Same default output column name
  `"nr_child_exp"` as `nr_children_experiment`, but a different number.**
- `.add_nr_children(...)` — internal; calls `nr_obs_sample` with an explicit
  distinct `new_child`.

### `R/LFQData.R`
- `.make_name_AinB(level_a, level_b, prefix = "nr_")` → name `nr_<B>_IN_<A>`.
- `.nr_B_in_A(data, level_a, level_b, merge = TRUE)` — **distinct count** of
  level_b per level_a (e.g. peptides per protein); merges back, `message()`s the
  added column name. Not observation-aware.
- `nr_B_in_A(pdata, config, merge = TRUE)` — config-driven wrapper using
  `hierarchy_keys_depth()` as A and the next hierarchy level as B. Used by
  `filter_proteins_by_peptide_count()`.

### `R/tidyMS_summarize_hierarchy.R`
- `hierarchy_counts(pdata, hierarchy_keys, isotope_label)` — `n_distinct` of
  each hierarchy level grouped by isotope label (experiment-wide, no per-sample
  split, no observed filter).
- `summarize_hierarchy(pdata, hierarchy_keys, isotope_label, hierarchy, factors)`
  — `n_distinct` of the other hierarchy levels grouped by a chosen level (e.g.
  peptides per protein); output columns named `<level>_n`.
- `HierarchyCountsSample` R6 + `hierarchy_counts_sample(...)` — `n_distinct` per
  isotope×sample, filtered by response not-NA **and** `nr_children >= threshold`.

## The actual problems

1. **At least four functions compute "peptides per protein"** with different
   conventions: `nr_B_in_A`, `summarize_hierarchy`, `nr_features_experiment`,
   and (observation-aware) `nr_children_experiment`. No single documented place
   says which to use when.
2. **Name collision**: `nr_children_experiment` and `nr_features_experiment`
   default their output column to the **same** `"nr_child_exp"` while computing
   semantically different quantities (observed-summed-max vs distinct count).
3. **Three naming schemes** for the same idea: `nr_<B>_IN_<A>` (`nr_B_in_A`),
   `<level>_n` (`summarize_hierarchy`), bare level name (`hierarchy_counts`),
   and `nr_child_exp` (the `*_experiment` pair).
4. **Misleading names**: `nr_obs_sample` sums a count column rather than
   counting observations; "children" / "features" / "B_in_A" all denote the
   same hierarchy relationship.
5. **Silent observation-awareness split**: some count distinct identifiers
   regardless of whether measured; others count only non-NA-response rows.
   Choosing the wrong one changes weights / DEqMS counts without any error.
6. **Scattered across three files** with no shared taxonomy or `@family` tag.

## Open questions (decide before refactoring)

- What are the orthogonal axes? Tentatively: **(a)** count-distinct vs
  sum-existing-count; **(b)** observation-aware vs identifier-only; **(c)**
  per-sample vs per-isotope vs experiment-wide; **(d)** reduce-across-samples
  (max/sum) or not.
- Can these collapse into one parameterized counter (e.g.
  `count_children(data, a, b, observed =, per_sample =, reduce =)`) plus thin
  named wrappers, without breaking the public/exported surface?
- Which of these are exported and depended on downstream (`prolfquapp`,
  `prophosqua`, `prolfquasaint`)? Grep before changing any signature — several
  are `@export`/`@keywords internal`.
- At minimum (low-risk, do first): give `nr_features_experiment` a distinct
  default output column name to remove the `"nr_child_exp"` collision, and add a
  cross-referencing `@family` doc block listing the whole counting family so the
  intended choice is discoverable.

## Constraints

- Do not change exported signatures without checking downstream callers and the
  `prolfqua` R6 class/API contract.
- Add a regression test pinning each counter's semantics before consolidating
  (the `*_experiment` pair already has `tests/testthat/test-tidyMS_aggregation.R`).
