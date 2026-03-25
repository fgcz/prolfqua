## Performance Code Review: prolfqua R Package — All Items Complete

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

### ~~5. `contrasts_linfct` per-protein loop + per-row `data.frame()` creation~~ DONE
Fixed: `compute_contrast_vectorized()` uses matrix multiplication (`linfct %*% coef`) with zero-filled coefficients and invalid-row masking for NA coefficients. Standard errors via `sqrt(diag(L %*% Sigma %*% t(L)))`. Original `compute_contrast()` untouched; hot-swap via `options(prolfqua.vectorize = TRUE)`.

---

### ~~6. `rlang::parse_expr()` + `mutate()` in a loop~~ DONE
Fixed: `linfct_matrix_contrasts_vectorized()` pre-parses all expressions with `lapply(parse_expr)` and evaluates in a single `dplyr::mutate(data, !!!parsed)` call. Falls back to per-expression loop on error. Original `linfct_matrix_contrasts()` untouched; hot-swap via `options(prolfqua.vectorize = TRUE)`.

---

### ~~7. `pivot_model_contrasts_2_Wide` uses `reduce(left_join)`~~ DONE
Fixed: replaced loop of N `pivot_wider` + `reduce(left_join)` with single `pivot_wider(values_from = columns, names_glue = "{.value}.{contrast}")`. Eliminates the `m_spread` helper, the loop, and the `purrr::reduce`/`dplyr::left_join` chain.

---

### ~~8. `apply()` instead of `rowSums` / partial sort~~ DONE
Fixed: `estimate_lod_global()` uses `rowSums(is.na(...))` instead of `apply()`. `function_lod_quantile()` keeps data as matrix (no data.frame conversion) and uses `sort(x, partial = n)[1:n]` for O(n) partial sort.

---

### ~~9. `get_LOD()` not cached~~ DONE
Fixed: added `private$.lod_cache` field to `MissingHelpers`. `get_LOD()` computes once on first call, returns cached value thereafter. Same null-check pattern already used for `self$stats`.

---

### ~~10. Filter-compute-rejoin pattern~~ DONE
Fixed: replaced filter → compute → select → `left_join` with a single `mutate` using `purrr::map2_*` that checks `has_model_fit` per row. Failed fits get `NA` directly from the else branch.

---

### All 10 items complete.
