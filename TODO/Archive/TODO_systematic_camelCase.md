# Plan: Systematic camelCase audit and migration strategy for prolfqua

## Context

The `AnalysisConfiguration` field migration (camelCase → snake_case) is complete. Now the user wants to address camelCase systematically across the **entire prolfqua package**. This plan is an analysis and strategy document.

COMMENT: Aim is to use snake_case variable names througout. Dataframe column names stay as is!

## Current State: Naming Convention Audit

### Overall ratio
- **~75% snake_case, ~25% camelCase** across the codebase
- R6 method names: 95% snake_case already
- Standalone function names: 90% snake_case already
- Function **parameters**: most problematic — many camelCase survivors
- Local variables inside function bodies: mixed
- Data column name defaults: camelCase (`"sampleName"`, `"isotopeLabel"`)



---

## Category 1: R6 Method Names (LOW count, HIGH impact)

### Deprecated aliases (callers exist — fix callers, then delete)
| Method | Snake_case equivalent | Callers |
|--------|----------------------|---------|
| `hierarchyKeys()` | `hierarchy_keys()` | `LFQData.R:385`, prolfquapp (2 files) |
| `hkeysDepth()` | `hierarchy_keys_depth()` | `LFQData.R:384` |

COMMENT: go for it.

### Factory methods using `get_CamelCase()` pattern (26 classes, ~6 methods)
These are the most visible camelCase in the API:
| Current | Proposed |
|---------|----------|
| `get_Transformer()` | `get_transformer()` |
| `get_Aggregator()` | `get_aggregator()` |
| `get_Plotter()` | `get_plotter()` |
| `get_Stats()` | `get_stats()` |
| `get_Summariser()` | `get_summariser()` |
| `get_Imputer()` | `get_imputer()` |

Called from: prolfqua internally, prolfquapp, prophosqua, prolfquasaint, prolfquabenchmark — **very high blast radius**.

COMMENT: keep as is these are Class names.

### Other camelCase R6 methods

| Method | Class | File |
|--------|-------|------|
| `subject_Id()` | `LFQData` | LFQData.R |
| `omit_NA()` | `LFQData` | LFQData.R |
| `NA_heatmap()` | `LFQDataPlotter` | LFQDataPlotter.R |
| `get_LOD()` | `MissingHelpers` | tidyMS_missingness_imputation.R |

COMMENT: rename - we anyway have a lot of breaking changes.

---

## Category 2: R6 Field Names (HIGH count, HIGH impact)

Fields accessed by downstream packages via `$fieldName`:

COMMENT: we have to minimize the access by field name. Access must be through methods because long term plan might be improving representation with LFQData. Please do not mix up data frame column names and variable names!


| Field | Classes using it | Proposed |
|-------|-----------------|----------|
| `modelName` | Contrasts, ContrastsPlotter, ContrastsTable, ContrastsMissing, ContrastsModerated, ContrastsFirth, ContrastsROPECA, Model, ModelFirth | `model_name` |
| `subject_Id` | Contrasts*, Model*, ContrastsPlotter | `subject_id` |
| `contrastDF` | ContrastsPlotter | `contrast_df` |
| `modelDF` | Model | `model_df` |
| `p.adjust` | Contrasts, ContrastsMissing, ContrastsModerated, ContrastsFirth, ContrastsROPECA, Model | Already uses `.` not camelCase — R convention for stats |
| `avg.abundance` | ContrastsPlotter | Same — R stats convention |
| `protein_annot` | Contrasts, ContrastsPlotter | Already snake_case |

**Note:** `p.adjust` and `avg.abundance` use `.` separator — this follows R stats conventions (like `p.adjust()` in base R). Probably leave as-is.

---

## Category 3: Exported Function Names (11 functions)

| Function | File | Proposed |
|----------|------|----------|
| `scriptCopyHelperVec` | vignetteHelpers.R | `script_copy_helper_vec` |
| `pivot_model_contrasts_2_Wide` | tidyMS_contrasts.R | `pivot_model_contrasts_to_wide` |
| `get_UniprotID_from_fasta_header` | utilities.R | `get_uniprot_id_from_fasta_header` |
| `isSingular_lm` | tidyMS_build_model.R | `is_singular_lm` |
| `plot_NA_heatmap` | tidyMS_plotting.R | OK — `NA` is an acronym |
| `sim_lfq_data_2Factor_config` | simulate_LFQ_data.R | `sim_lfq_data_2factor_config` |
| `remove_NA_rows` | utilities.R | OK — `NA` is an acronym |
| `nr_B_in_A` | LFQData.R | already snake + abbreviation |
| `list_to_AnalysisConfiguration` | AnalysisConfiguration.R | OK — class name is PascalCase |
| `aggregate_intensity_topN` | tidyMS_aggregation.R | OK or `aggregate_intensity_top_n` |
| `UpSet_interaction_missing_stats` | tidyMS_summarize_hierarchy.R | `upset_interaction_missing_stats` |


Comment: get_UniprotID_from_fasta_header where do we use it, I really would like to drop or move it to prolfquapp.
---

## Category 4: Exported Function Parameters (~19 functions)

Most impactful camelCase parameters:

| Parameter | Used in | Proposed |
|-----------|---------|----------|
| `colName` | `set_response()` | `col_name` |
| `workIntensity` | `make_reduced_hierarchy_config()` | `work_intensity` |
| `modelName` | `build_model_limma()`, `build_model_limpa()`, `merge_contrasts_results()`, ContrastsLimma init, ContrastsROPECA init, `plot_lrt_diagnostics()` | `model_name` |
| `preserveMean` | `robust_scale()`, `robscale()`, `robscale_subset()`, `scale_with_subset()`, `transform_work_intensity()` | `preserve_mean` |
| `sampleName` | `medpolish_estimate()`, `medpolish_estimate_df()` | `sample_name` |
| `factorDepth` | `omit_NA()` | `factor_depth` |
| `proteinName` | `plot_hierarchies_line()` | `protein_name` |
| `intesityNewName` | `transform_work_intensity()` | `intensity_new_name` (also fix typo) |
| `modelDF` | `compute_borrowed_variance()` | `model_df` |

---

## Category 5: Data Column Name Defaults

The config stores default **data column names** as strings:
```r
sample_name = "sampleName"       # config field is snake_case, but default column name is camelCase
isotope_label = "isotopeLabel"   # same pattern
```

These propagate into actual data frames via `setup_analysis()`. Changing these would break **all existing datasets** and serialized configs.

---

## Category 6: Local Variables Inside Function Bodies

Hundreds of camelCase local variables (`valueVars`, `annotationVars`, `idVars`, etc.). These are internal and don't affect the public API.

---

## Recommended Strategy: Phased Migration

### Phase 0: Quick wins — deprecated aliases (DO NOW)
- Fix callers of `hierarchyKeys()` / `hkeysDepth()` (4 call sites)
- Fix `AnalysisTableAnnotation$new()` in prolfquasaint (2 files)
- Delete the deprecated aliases from AnalysisConfiguration
- **Blast radius:** Small, self-contained

### Phase 1: R6 field names `modelName` / `subject_Id` / `contrastDF`
- These are used in output data frames and consumed downstream
- Add active bindings for backward compat? Or hard-switch?
- **Blast radius:** Large — touches all Contrasts/Model classes + downstream

COMMENT: data frame column names remain untouched, you have to be carefull, not to mix variable names and dataframe column names, unfortunately access looks same, $name. Shall we possibly make it clear by a naming convention, when do we operate on a dataframe and when on a class, e.g. prefix or suffix variables with data frames with df_ or _df?



### Phase 2: Factory method names `get_Transformer()` → `get_transformer()`
- These are the most visible API surface
- Add deprecated forwarding wrappers (like `hierarchyKeys` pattern)
- Migrate all internal + downstream callers
- Remove wrappers
- **Blast radius:** Very large — every analysis script uses these

COMMENT: do not touch see above.

### Phase 3: Exported function parameters
- `modelName` → `model_name`, `preserveMean` → `preserve_mean`, etc.
- Can't use active bindings for parameters — need to accept both old + new during transition, or hard-switch
- **Blast radius:** Medium — mostly internal usage

COMMENT: we go for hard switch.

### Phase 4: Exported function names
- ~6 functions to rename
- Add deprecated wrappers that call the new name + emit `.Deprecated()` warning
- **Blast radius:** Medium

COMMENt: we go for hard switch.

### Phase 5: Local variables (OPTIONAL / LOW PRIORITY)
- Purely cosmetic — no API impact
- Do opportunistically when touching a file for other reasons

COMMENT: yes local variables certainly should go to snake_case, I even would say we should tackle this first!

### DO NOT CHANGE: Data column name defaults
- `"sampleName"`, `"isotopeLabel"` as default column values — changing these breaks stored data and all existing workflows

---

## Decision Points for the User

1. **Hard-switch vs deprecation wrappers?** Wrappers are safer but add code. Given the ecosystem is small and you control all packages, a hard-switch + grep-fix-all-downstream may be faster.

2. **Which phases to do now?** Phase 0 is ready. Phases 1–4 could be done incrementally or as one batch.

3. **`get_Transformer()` naming:** Is the convention `get_transformer()` or just `transformer()`? Some R6 APIs drop the `get_` prefix entirely.

4. **`p.adjust` / `avg.abundance`:** Leave the dot-separated names? They follow R stats conventions.

5. **Scope:** Do you want to address all of Phases 0–4, or focus on specific categories?


Comment: Overall I would go for hard switch. Commit first all current changes. and then go hard. But also accompany it by encapsulating all class varibles -> We can even think going private. This also means, all methods with signature function(data, config), moving them ideally to function(LFQData) and then accessing the fields they need through methods.