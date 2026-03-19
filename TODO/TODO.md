# prolfqua — Open TODO Items

---

## 1. Limma: Expose Remaining lmFit Paths

From `Expose_additional_lmFit_paths.md`. Weights support is done; two features remain:

- **`block` / `duplicateCorrelation`** — for repeated measures / paired sample designs. Requires calling `limma::duplicateCorrelation()` first, then passing `block=` + `correlation=` to `lmFit()`. Needs `statmod` package (add to Suggests).
- **`method = "robust"`** — robust fitting via `MASS::rlm` (mrlm). Down-weights outlier observations at the fitting stage (separate from `robust = TRUE` in `eBayes()`). Pass `method` through to `lmFit()`.

Files: `R/ContrastsLimma.R`, `tests/testthat/test-ContrastsLimma.R`

---

## 2. API Surface Cleanup

**Audit completed (2026-03-19):** Full audit of all 113 functions with both `@export` and `@keywords internal`. Detailed results with exact caller locations across all 5 downstream packages in `TODO/API_surface_audit.md`.

**Findings:**
- **11 functions** called directly as `prolfqua::fn()` by downstream packages (Tier 1)
- **100 functions** internal (Tier 2) — includes 5 functions reached downstream only via R6 methods (e.g. `lfqdata$remove_small_intensities()`) — the R6 class is the public API, not the free function
- **2 S3 methods** (`sigma.rlm`, `df.residual.rlm`) — candidates for `@exportS3Method`

**Done (2026-03-19):** Removed `@keywords internal` from 9 Tier 1 functions. Two (`linfct_from_model`, `linfct_matrix_contrasts`) kept internal — only called from `inst/MyArticle/`. Added missing `@param` tags. R CMD check passes (0 errors, 0 new warnings).

---

## 3. Test Coverage Expansion

From `TODO_example_redundancy_reduction.md`. Three formal test files exist (`test-Model.R`, `test-ContrastsFacades.R`, `test-ContrastsLimma.R`). Still missing:

- `tests/testthat/test-LFQData.R` — LFQData + decorator classes (Transformer, Stats, Aggregator, Plotter, Summariser, Imputer)
- `tests/testthat/test-Contrasts.R` — All Contrasts* implementations (Wald, Moderated, ROPECA, Missing, Firth, Table)
- `tests/testthat/test-ContrastsPlotter.R` — ContrastsPlotter plot type verification

See the detailed plan previously in `TODO_example_redundancy_reduction.md` for per-file test specifications.

---

## 4. Example Runtime Reduction

From `TODO_example_redundancy_reduction.md`. R CMD check examples take ~70s. Opportunities:

- **15+ files** create identical `sim_lfq_data_peptide_config()` setup redundantly
- **5+ files** load `prolfqua_data('data_ionstar')$filtered()` (heavy)
- `LFQDataTransformer.R` loads `data_ionstar` twice in two separate example blocks
- `tidyMS_plotting.R` has 5 free function examples that duplicate `LFQDataPlotter` class examples
- Reduce `Nprot` from 100 to 20 in Contrasts examples
- Wrap heavy optional examples in `\donttest{}` (`LFQDataToSummarizedExperiment`, `Benchmark`)

Target: reduce from ~70s to ~50-55s.

---

## 5. Relocate Debugging Tips (Optional)

The deleted `SKILL.md` contained useful debugging tips (vignette variable shadowing in `globalenv()`, missing `stats::` namespace prefix, `Rplots.pdf` artifact). Consider moving these to:

- `CLAUDE.md` (so Claude Code has them in context), or
- A developer guide / contributing doc

---

## 6. Code Review: Remaining Minor Items

Open items from the 2026-02-19 code review not yet addressed:

- Dead `NA <- NA` assignment in `R/LFQDataImp.R`
- Empty `LFQDataImp` class body (`R/LFQDataImp.R:152-165`)
- `Contrasts_from_factors.Rmd` is a stub vignette (18 lines)
- Mixed `for()` loops and `map()` in `intensity_summary_by_hkeys()` (`R/tidyMS_aggregation.R:866-895`)
- Missing `@param` documentation on some utility functions in `R/utilities.R`

---

## Cross-Project Reference

The integration test infrastructure plan (previously `dazzling-greeting-nygaard.md`) belongs in `prolfqua_fml/integration_test/TODO/` — see `TODO_integration_test.md` in that directory.
