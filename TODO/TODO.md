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

## 4. R6 Encapsulation Follow-Up

The LFQData `data`/`config` encapsulation work is complete and the detailed
analysis note has been archived at
`TODO/Archive/TODO_LFQData_access_patterns.md`.

Still to consider as separate follow-up work:
- AnalysisConfiguration: `add_hierarchy()`, `add_factor()` methods for structured fields
- Contrasts/Model: make `model_df`, `model_name`, `subject_id` private

---

## Cross-Project Reference

The integration test infrastructure plan is in `prolfqua_fml/integration_test/TODO/TODO_integration_test.md`.
