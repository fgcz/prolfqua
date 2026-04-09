# LFQData Access Pattern Analysis

## Purpose

Audit every function that takes `(data, config)` and every decorator that accesses `self$lfq$data` / `self$lfq$config` directly. The goal: define a minimal LFQData API so that (a) standalone functions receive only what they need, and (b) decorators never touch raw fields — only methods that return data frames, vectors, or scalars.

Once the API is defined, item 4 (R6 Encapsulation) becomes mechanical: make `data` and `config` private, expose them through these methods.

---

## 1. Config Access Patterns (AnalysisConfiguration)

Every `(data, config)` function and every decorator accesses config. The accesses cluster into **5 roles**:

### Role A: Response column management

| Accessor | Returns | Used by |
|----------|---------|---------|
| `config$get_response()` | `character(1)` — current intensity column name | transform_work_intensity, apply_to_response_matrix, scale_with_subset, medpolish_estimate_dfconfig, rlm_estimate_dfconfig, aggregate_intensity_top_n, nr_obs_sample, summarize_stats, encode_bin_resp, remove_small_intensities, plot_intensity_distribution_*, rank_peptide_by_intensity |
| `config$set_response(name)` | void | LFQDataTransformer (after adding transformed column) |
| `config$pop_response()` | `character(1)` | LFQDataTransformer (robscale renames) |
| `config$is_response_transformed` | `logical(1)` | summarize_stats, plot_intensity_distribution_*, LFQDataStats, Aggregators |

**LFQData API mapping:** `response()` already exists. `is_transformed()` already exists. The set/pop are only used by Transformer which operates on a clone — keep as internal mutation on the clone's config.

### Role B: Hierarchy (measurement levels)

| Accessor | Returns | Used by |
|----------|---------|---------|
| `config$hierarchy_keys()` | `character(n)` — all hierarchy column names | complete_cases, upset_missing_stats, nr_obs_experiment, hierarchy_counts, summarize_hierarchy, summarize_stats |
| `config$relevant_hierarchy_keys()` | `character(m)` — hierarchy columns at current depth | medpolish_estimate_dfconfig, rlm_estimate_dfconfig, estimate_intensity, plot_estimate, aggregate_intensity_top_n, nr_obs_sample, nr_obs_experiment, nr_B_in_A, sample_subset, separate_hierarchy, plot_hierarchies_line_df, LFQDataPlotter, AggregateLimpa |
| `config$hierarchy` | named list (full hierarchy definition) | separate_hierarchy |
| `config$hierarchy_depth` | `integer(1)` | (set in examples, read by hierarchy_keys_depth) |
| `config$sep` | `character(1)` — hierarchy separator | separate_hierarchy |

**LFQData API:** `subject_id()` already wraps `relevant_hierarchy_keys()` (rename from `subject_Id()`). Keep `relevant_hierarchy_keys()` name — no obviously better alternative.

### Role C: Sample / file identity

| Accessor | Returns | Used by |
|----------|---------|---------|
| `config$sample_name` | `character(1)` — sample column name | complete_cases, .reestablish_condition, plot_pca, plot_heatmap, plot_heatmap_cor, plot_raster, plot_na_heatmap, plot_intensity_distribution_*, plot_sample_correlation, nr_obs_sample, aggregate_intensity_top_n, response_matrix_as_tibble |
| `config$file_name` | `character(1)` — file column name | complete_cases, .reestablish_condition, aggregate_intensity_top_n, nr_obs_sample |
| `config$isotope_label` | `character(1)` or `NULL` | complete_cases, .reestablish_condition, summarize_stats, hierarchy_counts, summarize_hierarchy, response_matrix_as_tibble |

**LFQData API:** Need `sample_name()`, `file_name()` accessors. `isotope_label` — expose but document as legacy.

### Role D: Experimental factors

| Accessor | Returns | Used by |
|----------|---------|---------|
| `config$factor_keys()` | `character(n)` — all factor column names | complete_cases, .reestablish_condition, plot_sample_correlation, plot_heatmap_cor, plot_heatmap, plot_raster, plot_na_heatmap, plot_pca, aggregate_intensity_top_n |
| `config$relevant_factor_keys()` | `character(m)` — factor columns at current depth | summarize_stats, plot_stat_violin, plot_stat_violin_median, plot_stdv_vs_mean, LFQDataStats |
| `config$factor_depth` | `integer(1)` | summarize_stats_factors, LFQDataStats (modifies!) |

**LFQData API naming decision:**
- `config$factor_keys()` → expose as `factor_keys()` on LFQData
- `config$relevant_factor_keys()` → expose as `relevant_factor_keys()` on LFQData

### Role E: Miscellaneous metadata

| Accessor | Returns | Used by |
|----------|---------|---------|
| `config$nr_children` | `character(1)` — name of nr_children column | .add_nr_children, nr_obs_sample, nr_obs_experiment, aggregate_intensity_top_n |
| `config$work_intensity` | `character(1)` | Aggregators (warning messages only) |
| `config$bin_resp` | computed field | encode_bin_resp |
| `config$ident_q_value` | `character(1)` | aggregate_intensity_top_n |
| `config$complete_cases()` | function | upset_missing_stats |

---

## 2. Data Access Patterns

Functions access `data` in these ways:

### Read-only (majority)
Most functions receive `data` as a data frame, select/filter/group/summarize using column names from config, and return a new data frame. They never mutate the input.

Examples: all plotting functions, summarize_stats, hierarchy_counts, complete_cases, rank_peptide_by_intensity

### Read + transform → new data frame
Functions that produce a modified copy: transform_work_intensity, apply_to_response_matrix, scale_with_subset, estimate_intensity, aggregate_intensity_top_n

**Transformer handling:** These functions return a new data frame. The Transformer decorator currently assigns the result back to `self$lfq$data`. Under option 3 (return new LFQData — decided), each transform step creates a new LFQData instead of mutating the clone. No setters needed on LFQData — Transformer constructs new instances internally.

### Read + join → augmented data frame
Functions that add columns: encode_bin_resp, .add_nr_children, nr_obs_sample

**Suggestion:** Same pattern as transform — these return augmented data frames. Callers (LFQData methods or decorators) can construct a new LFQData from the result if needed. No special API required.

---

## 3. Decorator Access Patterns

### Current pattern (problematic)
All decorators (except Summariser) clone the LFQData, then reach through to `self$lfq$data` and `self$lfq$config` to call standalone functions:

```r
# Typical decorator method
plot_density = function() {
  plot_intensity_distribution_density(self$lfq$data, self$lfq$config)
}
```

This means `data` and `config` can never be made private.

### Decorator-specific access summary

| Decorator | Reads data | Writes data | Reads config | Writes config | Uses LFQData methods |
|-----------|-----------|-------------|-------------|---------------|---------------------|
| **Transformer** | yes | **yes** | yes | **yes** (set/pop_response) | is_transformed(), clone() |
| **Plotter** | yes | constructor only (NA removal) | yes | no | clone() |
| **Stats** | yes | no | yes | **yes** (factor_depth) | is_transformed(), clone() |
| **Summariser** | yes | no | yes | no | hierarchy_counts() (LFQData method) |
| **Aggregators** | yes | no | yes | no | is_transformed(), clone(), to_wide() |

**Key insight:** Only Transformer and Stats mutate the internal state. They operate on *clones*, so the mutation is safe.

**Transformer decision:** Option 3 — return new LFQData from each transform step instead of mutating the clone. This is cleaner and avoids needing write access to LFQData fields entirely.

---

## 4. Proposed LFQData API

### Tier 1: Column name accessors (return character vectors)

```r
# Already exist (with renames)
response()              # character(1) — intensity column name
subject_id()            # character(m) — hierarchy columns at current depth (rename from subject_Id)

# New accessors
hierarchy_keys()            # character(n) — all hierarchy column names
relevant_hierarchy_keys()      # character(m) — hierarchy columns at current depth (same as subject_id)
factor_keys()               # character(n) — all factor column names (wraps config$factor_keys)
relevant_factor_keys()      # character(m) — factor columns at current depth (wraps config$factor_keys_depth)
sample_name()           # character(1) — sample column name
file_name()             # character(1) — file/raw-file column name
nr_children_col()       # character(1) — name of the nr_children column
```

### Tier 2: Data accessors (return data frames)

```r
# Already exist
to_wide(...)            # list(data, annotation, rowdata, config) — wide format

# New accessor
get_data()              # data.frame — the tidy intensity data
```

Note: No `get_config()`. Config is not exposed. All config access goes through the API methods above.

### Tier 3: State queries (return scalars)

```r
# Already exist
is_transformed()        # logical(1) — is intensity log-transformed?
```

No `hierarchy_depth()` or `factor_depth()` scalars — derivable from `length(relevant_hierarchy_keys())` etc.

### Tier 4: Derived data (return data frames)

```r
# Already exist
hierarchy_counts()      # data.frame — count of items at each hierarchy level
factors()               # data.frame — unique factor combinations per sample (no collision: factor_keys() returns column names)
```

---

## 5. Per-Function Config Access Detail

### Group 1: Plotting functions — tidyMS_plotting.R

These are called exclusively from LFQDataPlotter. Config is NOT to be exposed — each function should receive only what it needs.

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `plot_intensity_distribution_violin` | `sample_name`, `get_response()`, `is_response_transformed` | `sample_col`, `response_col`, `is_transformed` |
| `plot_intensity_distribution_density` | `get_response()`, `sample_name`, `is_response_transformed` | `response_col`, `sample_col`, `is_transformed` |
| `plot_sample_correlation` | indirect via `tidy_to_wide_config()` | needs `tidy_to_wide_config` refactored first |
| `plot_heatmap_cor` | `factor_keys()`, `sample_name` + `tidy_to_wide_config()` | `factor_cols`, `sample_col` + wide matrix |
| `plot_heatmap` | `factor_keys()`, `sample_name` + `tidy_to_wide_config()` | `factor_cols`, `sample_col` + wide matrix |
| `plot_raster` | `factor_keys()`, `sample_name` + `tidy_to_wide_config()` | `factor_cols`, `sample_col` + wide matrix |
| `plot_na_heatmap` | `factor_keys()`, `sample_name` + `tidy_to_wide_config()` | `factor_cols`, `sample_col` + wide matrix |
| `plot_pca` | `sample_name`, `factor_keys()[1]`, `factor_keys()[2]` + `tidy_to_wide_config()` | `sample_col`, `color_col`, `shape_col` + wide matrix |

**Pattern:** The heatmap/PCA functions all call `tidy_to_wide_config(data, config, as.matrix = TRUE)` first, then use `factor_keys()` and `sample_name` for annotation. Refactoring `tidy_to_wide_config` is the key dependency.

**Proposed refactoring:** The Plotter decorator calls `self$lfq$to_wide(as.matrix = TRUE)` to get the wide matrix + annotation, then passes the pre-computed pieces to the plotting functions. The plotting functions no longer need config at all — they get `(matrix, annotation, ...)`.

### Group 2: Plotting functions — tidyMS_stats.R

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `plot_stat_density` | `factor_keys()[1]` | `color_col` — single string |
| `plot_stat_density_median` | `factor_keys()[1]` | `facet_col` — single string |
| `plot_stat_violin` | `relevant_factor_keys()` | `group_cols` — character vector for `tidyr::unite()` |
| `plot_stat_violin_median` | `factor_keys()[1]` | `x_col` — single string |
| `plot_stdv_vs_mean` | `relevant_factor_keys()` | `facet_cols` — character vector |

**Pattern:** These only need 1-2 column names from config. Trivial to refactor.

### Group 3: Plotting functions — tidyMS_aggregation.R

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `plot_hierarchies_line_df` | `relevant_hierarchy_keys()` + passes config to `plot_hierarchies_line` | `hierarchy_cols` + all of plot_hierarchies_line's needs |
| `plot_hierarchies_line` | `hierarchy_keys(TRUE)`, `sample_name`, `get_response()`, `relevant_factor_keys()`, `isotope_label`, `is_response_transformed` | 6 individual args |
| `plot_estimate` | `relevant_hierarchy_keys()` (×2), `get_response()` on config_reduced | `hierarchy_cols`, `response_col` |

### Group 4: Missingness functions

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `missigness_histogram` | `relevant_factor_keys()`, `is_response_transformed`, `isotope_label`, `get_response()` | 4 individual args |
| `upset_missing_stats` | `get_response()`, `hierarchy_keys()`, `sample_name` | 3 individual args |

### Group 5: Stats functions — tidyMS_stats.R

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `summarize_stats` | `get_response()`, `hierarchy_keys()`, `isotope_label`, `relevant_factor_keys()`, `is_response_transformed`, `factor_keys()` | 6 fields — candidate for `(lfqdata)` |
| `summarize_stats_factors` | `factor_depth`, `relevant_factor_keys()` | 2 fields |
| `summarize_stats_all` | delegates to `summarize_stats` | same as above |

### Group 6: Data setup / hierarchy — tidyMS_data_setup.R, tidyMS_summarize_hierarchy.R

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `complete_cases` | `file_name`, `sample_name`, `factor_keys()`, `isotope_label`, `hierarchy_keys()` | 5 fields — candidate for `(lfqdata)` |
| `sample_subset` | `relevant_hierarchy_keys()` | 1 field |
| `separate_hierarchy` | `relevant_hierarchy_keys()`, `hierarchy` (named list), `sep` | 3 fields, needs raw hierarchy definition |
| `hierarchy_counts` | `hierarchy_keys()`, `isotope_label` | 2 fields |
| `summarize_hierarchy` | `isotope_label`, `hierarchy_keys()` | 2 fields |

### Group 7: Aggregation functions — tidyMS_aggregation.R

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `medpolish_estimate_dfconfig` | `hierarchy_keys()`, `relevant_hierarchy_keys()`, `get_response()`, `sample_name` | 4 fields |
| `rlm_estimate_dfconfig` | same as medpolish | 4 fields |
| `estimate_intensity` | `relevant_hierarchy_keys()` | 1 field |
| `aggregate_intensity_top_n` | `relevant_hierarchy_keys()`, `sample_name`, `file_name`, `isotope_label`, `factor_keys()`, `get_response()`, `ident_q_value` | 7 fields — heaviest user, candidate for `(lfqdata)` |
| `.reestablish_condition` | `sample_name`, `factor_keys()`, `file_name`, `isotope_label` | 4 fields |
| `.add_nr_children` | `relevant_hierarchy_keys()`, `file_name`, `nr_children` | 3 fields |
| `nr_obs_sample` | `get_response()`, `relevant_hierarchy_keys()`, `file_name`, `nr_children` | 4 fields |
| `nr_obs_experiment` | `hierarchy_keys()`, `relevant_hierarchy_keys()`, `nr_children` | 3 fields |
| `rank_peptide_by_intensity` | `get_response()` | 1 field |

### Group 8: Transformer helpers — LFQDataTransformer.R

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `transform_work_intensity` | `get_response()`, `set_response()`, `is_response_transformed` | response management (Transformer-internal) |
| `response_matrix_as_tibble` | `sample_name`, `hierarchy_keys()`, `isotope_label`, `set_response()` | 4 fields |
| `get_robscales` | indirect via `tidy_to_wide_config()` | wide conversion |
| `apply_to_response_matrix` | `get_response()` | 1 field |
| `scale_with_subset` | `get_response()` | 1 field |

### Group 9: LFQData.R standalone functions

| Function | Config accesses | What it actually needs |
|----------|----------------|----------------------|
| `remove_small_intensities` | `get_response()` | 1 field |
| `nr_B_in_A` | `relevant_hierarchy_keys()`, `hierarchy_keys()` | 2 fields |

---

## 6. Refactoring Strategy

### Phase 1+3: Add API methods & update decorators — DONE (2026-04-08)

Added 8 API methods to LFQData (`hierarchy_keys()`, `relevant_hierarchy_keys()`, `factor_keys()`, `relevant_factor_keys()`, `sample_name()`, `file_name()`, `nr_children_col()`, `get_data()`).

Updated decorators to use API methods where they accessed config fields directly:
- **LFQDataPlotter** — `config$sample_name` → `self$lfq$sample_name()`, hierarchy depth logic → `self$lfq$relevant_hierarchy_keys()`/`hierarchy_keys()`, `pairs_smooth` → `self$lfq$get_data()`/`sample_name()`
- **LFQDataStats** — `config$hierarchy_keys()` → `self$lfq$hierarchy_keys()`, `config$factor_keys_depth()` → `self$lfq$relevant_factor_keys()`
- **LFQDataAggregator** — `.check_aggregatable` → `lfq$hierarchy_keys()`/`relevant_hierarchy_keys()`, AggregateLimpa → `self$lfq$relevant_hierarchy_keys()`
- **LFQDataSummariser** — no changes (all pass-throughs to standalone functions)
- **LFQDataTransformer** — skipped (has writes, deferred to Phase 2)

Pass-throughs `func(self$lfq$data, self$lfq$config)` left unchanged — standalone function signatures haven't changed yet.

### Phase 2: Refactor Transformer to return new LFQData — NEXT

Each transform method (`log2()`, `robscale()`, `vsn()`, etc.) returns a new LFQData instance instead of mutating the clone. The `$lfq` field becomes the result of the last transformation. No setters needed on LFQData.

### Phase 4: Refactor standalone function signatures

For functions with 1-3 config fields: replace `config` with individual named arguments.
For functions with 4+ fields (complete_cases, aggregate_intensity_top_n, summarize_stats): accept `lfqdata` and call API methods internally.

### Phase 5: Make `data` and `config` private

Use active bindings that warn on direct access (deprecation period), then error.

---

## 7. Pending Rename: `aggregate_intensity_topN` → `aggregate_intensity_top_n`

This function is still camelCase. Rename to snake_case as part of this work.

---

## 8. Key Dependency: `tidy_to_wide_config()`

Multiple plotting functions (heatmaps, PCA, correlation) call `tidy_to_wide_config(data, config, ...)` which is the main config-heavy conversion. Refactoring this function — or having LFQData's `to_wide()` method serve the same purpose — is a prerequisite for removing config from plotting functions.

**Proposal:** The Plotter decorator calls `self$lfq$to_wide(as.matrix = TRUE)` once, then passes the resulting matrix + annotation to plotting functions. This eliminates the `tidy_to_wide_config` dependency in individual plot functions.
