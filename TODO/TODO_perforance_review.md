## Performance Code Review: prolfqua R Package

Here are the issues ranked by impact, from highest to lowest:

---

### ~~1. O(n²) copy-on-modify in `impute_refit_singular`~~ DONE
Fixed: results collected into list first, then bulk-assigned with `[idx] <-`.

---

### ~~2. `complete_cases()` called N+1 times in `summarize_stats_factors`~~ DONE
Fixed: added `.completed` parameter to `summarize_stats`/`summarize_stats_all`; `summarize_stats_factors` and `LFQDataStats` now call `complete_cases` once and pass `.completed = TRUE` downstream.

---

### ~~3. Sequential `mutate()` chains creating unnecessary copies~~ DONE
Fixed: combined into single `mutate()` calls in both `tidyMS_build_model.R` (5 mutates → 1) and `tidyMS_moderation.R` (3 mutates → 1).

---

### ~~4. Deep-cloning in all decorator constructors~~ REVISED — mostly false alarm
**Original claim was wrong.** Critical review showed all 4 cloning decorators are justified:
- `LFQDataTransformer`, `LFQDataAggregator`, `LFQDataImp` — mutate data
- `LFQDataStats` — mutates `config$factorDepth` (line 97-98), clone needed
- `LFQDataSummariser` — already correctly avoids cloning

**Bug found and fixed:** `LFQDataPlotter` did NOT clone but DID mutate the original via `na.omit()` (line 53), silently dropping NA rows from the caller's data. Fixed by adding `$clone(deep = TRUE)`.

---

### 5. `contrasts_linfct` per-protein loop + per-row `data.frame()` creation
**`tidyMS_contrasts.R` ~lines 348-375, 473-529**

`compute_contrast` creates a `data.frame()` per row of `linfct` in a loop, then `bind_rows()`. This is called once per protein. Fix: pre-allocate result vectors and construct a single data.frame, or vectorize the matrix multiplication across all linfct rows.

---

### 6. `rlang::parse_expr()` + `mutate()` in a loop
**`tidyMS_contrasts.R` ~lines 145-198**

```r
for (i in seq_along(contrasts)) {
  data <- dplyr::mutate(data, !!names(contrasts)[i] := !!rlang::parse_expr(contrasts[i]))
}
```
Each iteration copies the tibble and parses an expression. Fix: build all contrasts in a single `mutate()` using `!!!` splicing.

---

### 7. `pivot_model_contrasts_2_Wide` uses `reduce(left_join)`
**`tidyMS_contrasts.R` ~lines 450-454**

Loops over columns, does a separate `pivot_wider` per column, then chains N-1 `left_join` calls. Fix: single `pivot_wider(values_from = columns, names_glue = "{.value}.{contrast}")`.

---

### 8. `apply()` instead of `rowSums` / partial sort
**`LFQDataImp.R` ~lines 14-17:**
```r
apply(data_matrix, 1, function(x) { sum(is.na(x)) })
```
Fix: `rowSums(is.na(data_matrix))`

**`LFQDataImp.R` ~lines 37-66:** Converts matrix to data.frame, then `sort()` per column.
Fix: keep as matrix, use `sort(x, partial = n)[1:n]` for O(n) partial sort.

---

### 9. `get_LOD()` not cached
**`tidyMS_missingness_imputation.R` ~lines 68-73**

Called multiple times (from `impute_weighted_lod`, `impute_lod`, `get_contrasts`), each time filtering and summarizing `get_stats()`. Fix: cache the LOD scalar in a private field.

---

### 10. Filter-compute-rejoin pattern
**`tidyMS_build_model.R` ~lines 595-627**

Filters to models with fits, computes columns, then `left_join` back to the full table. Fix: use conditional `mutate` with `if_else`/`case_when` on the full table, avoiding the filter + rejoin round-trip.

---

### Priority for remaining items:
Issues **5** and **6** run in the per-protein hot path and are the most impactful remaining items. Issues **7**, **8**, **9**, **10** are moderate improvements.