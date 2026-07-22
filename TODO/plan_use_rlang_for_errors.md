# Plan: Typed errors via `rlang` (scoped)

## Goal

Replace raw `stop()` at prolfqua's public boundaries with **typed `rlang` conditions**, so failures
carry machine-readable classes that tests and downstream packages can catch. Fill the one real
validation gap: the `LFQData` mutators.

**Scope.** A targeted swap of a few boundary `stop()` calls to typed conditions, plus validation on
the `LFQData` mutators — using only `rlang`, which is already a dependency.

## Why this scope

- `rlang` is **already in `Imports`** (currently used only for `sym`/`parse_expr`). Adopting
  `rlang::abort()` costs zero new dependencies.
- The one benefit that matters is **error classes** you can assert on — not nicer message text.
- The genuine soft spot is `LFQData$set_data()` / `set_config_value()` ([R/LFQData.R:103,117](../R/LFQData.R#L103))
  mutating internal state with no validation.

## Steps

### 1. Internal helper layer (new file `R/conditions.R`)

Small, internal (not exported). One base class `prolfqua_error` with subclasses.

```r
abort_bad_argument <- function(arg, must, not = NULL, parent = NULL) {
  msg <- c("Invalid argument.", x = sprintf("`%s` must %s.", arg, must))
  if (!is.null(not)) msg <- c(msg, i = sprintf("Actual value: %s.", not))
  rlang::abort(msg, class = c("prolfqua_error_bad_argument", "prolfqua_error"), parent = parent)
}

abort_missing_columns <- function(cols, data_nm = "data") {
  rlang::abort(
    c("Required columns are missing.",
      x = sprintf("Missing from `%s`: %s.", data_nm, paste(cols, collapse = ", "))),
    class = c("prolfqua_error_missing_column", "prolfqua_error_bad_argument", "prolfqua_error")
  )
}

abort_invalid_config <- function(msg) {
  rlang::abort(c("Invalid configuration.", x = msg),
    class = c("prolfqua_error_invalid_configuration", "prolfqua_error"))
}
```

### 2. Convert boundary `stop()` calls (behaviour-preserving)

Only these public entry points — same failure conditions, same triggers, just typed. Internal
helper `stop()` calls stay as they are:

- `setup_analysis()` — [R/tidyMS_data_setup.R:60,63,83,178](../R/tidyMS_data_setup.R#L60)
  (missing `file_name`, missing file-name column, no factors, duplicate hierarchy-key/sample).
  Keep the `debug = TRUE` diagnostic path and the existing `warning()` calls unchanged.
- `build_contrast_analysis()` — the `method` check. Keep the existing registry lookup
  (`lookup_facade()`); do not invent `list_facades()`.

Leave the other ~50 `stop()` calls in internal helpers alone — they are not public boundaries.

### 3. Validate `LFQData` mutators (the real gap)

Add a private `validate_state()` to `LFQData` and call it from `set_data()` and `set_config_value()`:
confirm the new data still carries the columns implied by the current config
(`file_name`, sample name, factor + hierarchy keys). On failure, `abort_bad_argument()` /
`abort_missing_columns()`. Reconcile column-name accessors with the **actual**
`AnalysisConfiguration` API before writing.

### 4. Tests

Add to the existing suites (`test-tidyconfig_functions.R`, `test-LFQData.R`,
`test-ContrastsFacades.R`). Assert on **class**, not message:

```r
expect_error(setup_analysis(bad, cfg), class = "prolfqua_error_missing_column")
expect_error(lfq$set_data(bad_tbl), class = "prolfqua_error_bad_argument")
expect_error(build_contrast_analysis(lfq, "~ group_", list(), method = "nope"),
             class = "prolfqua_error_bad_argument")
```

### 5. Docs / changelog

- `make document` (no NAMESPACE edits — helpers are internal).
- Add a `NEWS.md` bullet: user-visible errors from `setup_analysis()`, `build_contrast_analysis()`,
  and `LFQData` mutators now carry `prolfqua_error*` classes.

## Verification

`make test`, then `make check-fast`.

## Risk

Low. Behaviour is preserved — same errors fire on the same inputs, only their class changes.
Downstream code that matched on message text is unaffected (messages are kept); code can now
additionally catch by class. The only new *behaviour* is validation on `LFQData` mutators, which
previously accepted invalid state silently.
