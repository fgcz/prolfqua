# prolfqua — Open TODO Items

---

## 1. Limma: Expose Remaining lmFit Paths

Weights support is done; two features remain:

- **`block` / `duplicateCorrelation`** — for repeated measures / paired sample designs. Requires calling `limma::duplicateCorrelation()` first, then passing `block=` + `correlation=` to `lmFit()`. Needs `statmod` package (add to Suggests).
- **`method = "robust"`** — robust fitting via `MASS::rlm` (mrlm). Down-weights outlier observations at the fitting stage (separate from `robust = TRUE` in `eBayes()`). Pass `method` through to `lmFit()`.

Files: `R/ContrastsLimma.R`, `tests/testthat/test-ContrastsLimma.R`

---

## 2. Test Coverage Expansion

All three test files now exist (`test-LFQData.R`, `test-Contrasts.R`, `test-ContrastsPlotter.R`). May still benefit from additional coverage — review current assertions and add edge cases.

---

## 3. Code Review: Remaining Minor Items

Open items from the 2026-02-19 code review:

- `Contrasts_from_factors.Rmd` — stub vignette, 151 lines but may need completion or removal
- Mixed `for()` loops and `map()` in `intensity_summary_by_hkeys()` (`R/tidyMS_aggregation.R`)
- Missing `@param` documentation on some utility functions in `R/utilities.R`

*Resolved:* `LFQDataImp` class removed, `NA <- NA` gone.


---

## 5. Refactor `(data, config)` → `(lfqdata)` Signatures (Deferred)

COMMENT: I want a change to the item. I want you to go through all the functions with the signature (data, config) and analyse the access patterns to LFQData, the function then must only get what they need. not the entire LFQData object. basically expose to them only the essentials. This relates to item 4 which I moved behind item 5. Because once we know what are the access patterns, we can encapsulate, and design an interface, servicing the access patterns. 
Write down you analysis, with interface suggestions for LFQData in a file TODO/TODO_LFQData_access_patterns.md.
All decorators, that is, LFQPlotter, LFQAggregator, etc, must access the LFQData object throught the API. The API return types are data frames, arrays, strings, or lists.



~40 functions take separate `data` and `config` parameters. Many are only called from LFQData class methods that pass `self$data, self$config`. Refactor to accept an `lfqdata` object directly. Large effort, moderate benefit.

---

## 4. R6 Encapsulation (Deferred)

Add accessor methods and make fields private:

- AnalysisConfiguration: `add_hierarchy()`, `add_factor()` methods for structured fields
- LFQData: `get_data()` / `get_config()` methods, then make `data`/`config` private
- Contrasts/Model: make `model_df`, `model_name`, `subject_id` private (external code should use `get_contrasts()`, `get_coefficients()`, etc.)

---

## Cross-Project Reference

The integration test infrastructure plan is in `prolfqua_fml/integration_test/TODO/TODO_integration_test.md`.
