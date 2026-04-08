# Plan: Rename AnalysisConfiguration fields from camelCase to snake_case

## Context

AnalysisConfiguration is the central config class in prolfqua. Most fields are already snake_case, but 9 fields remained camelCase or mixed-case. Goal: normalize all field names to snake_case, with backward-compatible active binding aliases so downstream packages (prolfquapp, prophosqua, etc.) can migrate incrementally.

## Fields renamed

| Old | New | Default value (unchanged) |
|---|---|---|
| `fileName` | `file_name` | `NULL` |
| `sampleName` | `sample_name` | `"sampleName"` |
| `normValue` | `norm_value` | `NULL` |
| `isotopeLabel` | `isotope_label` | `"isotopeLabel"` |
| `workIntensity` | `work_intensity` | `NULL` |
| `factorDepth` | `factor_depth` | `1` |
| `hierarchyDepth` | `hierarchy_depth` | `1` |
| `ident_qValue` | `ident_q_value` | `"qValue"` |
| `ident_Score` | `ident_score` | `character()` |

**Important**: Default string values like `"sampleName"` and `"isotopeLabel"` are DATA column names, not config field names -- they stay unchanged.

## ALL STEPS COMPLETED

### ~~1. Rename fields + add active binding aliases~~ DONE

Renamed all 9 public fields to snake_case. Added 9 active bindings (no deprecation warnings) so old camelCase names still work for downstream compat. Updated all internal method references (`set_response`, `get_response`, `pop_response`, `factor_keys_depth`, `hierarchy_keys_depth`, `id_required`, `id_vars`, `value_vars`, `annotation_vars`).

### ~~2. Simplify the constructor~~ DONE

Dropped `analysisTableAnnotation` / `analysisParameter` parameters entirely. Constructor is now `initialize = function() {}` with no arguments.

### ~~3. Drop `table` and `parameter` active binding aliases~~ DONE

Removed the deprecated `table` and `parameter` active bindings. Fixed 2 internal callers.

### ~~4. Fix R6_extract_values() -- serialization~~ DONE

Added filtering via `.__enclos_env__$.__active__` to exclude active bindings from serialization output. Prevents duplicate keys in YAML/RDS.

### ~~5. Update make_reduced_hierarchy_config()~~ DONE

Renamed `workIntensity` -> `work_intensity` in function body.

### ~~6. Migrate prolfqua R/ files~~ DONE

All ~242 references across R/ files migrated. Data column names (`data$isotopeLabel`, `data$sampleName`, etc.) left unchanged.

### ~~7. Migrate tests, vignettes, data-raw, README~~ DONE

Updated references in tests/testthat/, vignettes/, data-raw/fix_deprecated_config.R, README.md, inst/issue71/.

### ~~8. Update CLAUDE.md field names~~ DONE

Updated architecture docs to use snake_case field names.

## What was NOT changed

- **Data column names** -- `"sampleName"`, `"isotopeLabel"`, `"qValue"` remain as column names in data
- **Downstream packages** -- prolfquapp, prophosqua, etc. keep working via active binding aliases; migrate separately
- **Method names** -- `set_response()`, `get_response()`, `hierarchy_keys()`, etc. already snake_case
- **Already-deprecated method aliases** -- `hierarchyKeys()` and `hkeysDepth()` stay as-is

## Downstream impact

- **Active binding aliases**: `config$fileName`, `config$sampleName`, etc. still work (read + write) via active bindings — no breakage for downstream packages
- **Constructor change**: `$new(atable)` pattern broken — ~17 calls across prolfquapp/prolfquappPTMreaders/prolfquasaint. Fix: drop the copy line, build config directly.
- **`config$table$X` / `config$parameter$X`**: now errors. Fix: `config$table$X` → `config$X`

All 618 tests pass (0 failures, 42 pre-existing warnings).
