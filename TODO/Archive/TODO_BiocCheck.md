# TODO: BiocCheck

Generated from `make check-bioc` with BiocCheck 1.46.3.

Current result:

```text
1 ERROR | 0 WARNINGS | 4 NOTES
```

## Error

- [ ] Ask the maintainer to add `prolfqua` to Watched Tags in the Bioconductor Support Site profile.
  - BiocCheck now recognizes the maintainer as registered at the Support Site.
  - Latest BiocCheck error: `Add package to Watched Tags in your Support Site profile`.
  - Visit: <https://support.bioconductor.org/accounts/edit/profile>

## Package Metadata Notes

- [ ] Confirm the maintainer is subscribed to the Bioc-Devel mailing list.
  - BiocCheck cannot determine this automatically because it requires admin credentials.
  - Subscribe/check here: <https://stat.ethz.ch/mailman/listinfo/bioc-devel>

## Coding Practice Notes

All coding-practice notes discussed above have been addressed. BiocCheck now reports `OK` for coding practice and parsed
R code checks.

## Documentation and Formatting Notes

- [x] Add runnable examples to exported-object documentation where feasible.
  - Addressed for the objects BiocCheck listed:
    `AnovaExtractor.Rd`, `build_model_limpa.Rd`, `FACADE_REGISTRY.Rd`, `HierarchyCountsSample.Rd`, `is_singular_lm.Rd`,
    `merge_contrasts_results.Rd`, `moderated_p_deqms_long.Rd`, `moderated_p_deqms.Rd`, `moderated_p_limma.Rd`,
    `panel_cor.Rd`, `pivot_model_contrasts_to_wide.Rd`, `plot_hierarchies_add_quantline.Rd`, `prolfqua_data.Rd`,
    `R6_extract_values.Rd`, `robust_scale.Rd`, `sample_subset.Rd`, `script_copy_helper_vec.Rd`, `strategy_limpa.Rd`,
    `StrategyLimma.Rd`, `StrategyLimpa.Rd`, `table_facade.Rd`, and `tidy_to_wide.Rd`.
  - Latest BiocCheck no longer reports missing examples.
- [ ] Review long functions; BiocCheck reports 37 functions over 50 lines.
  - Addressed: `fitFDistRobustly_LG()` in `R/squeezeVarRob.R` was covered by regression tests and split into helpers.
  - Direct golden-output tests now cover `fitFDist_LG()` plain, covariate-trend, missing-input, and validation paths, so
    it has a harness before any refactor.
  - Addressed: `fitFDist_LG()` in `R/squeezeVarRob.R` was split into internal helpers and is now 46 lines by
    `BiocCheck:::getFunctionLengths`; the remaining over-threshold function in that file is `trigammaInverse()` at 60
    lines.
- [ ] Decide whether to adopt Bioconductor formatting conventions over the package's local convention.
  - BiocCheck reports 1002 lines over 80 characters.
  - BiocCheck reports 5407 lines whose indentation is not a multiple of 4 spaces.
  - The tab note has been addressed; BiocCheck no longer reports tab-indented lines.
  - Do not bulk-reformat these remaining line-length/indentation notes without an explicit convention decision because the
    workspace AGENTS.md states 120-character lines and 2-space indentation.
