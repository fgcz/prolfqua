# TODO: BiocCheck

Generated from `make check-bioc` with BiocCheck 1.46.3.

Current result:

```text
1 ERROR | 3 WARNINGS | 21 NOTES
```

## Error

- [ ] Register or verify the maintainer email on the Bioconductor Support Site.
  - BiocCheck error: `Unable to find your email in the Support Site: HTTP 404 Not Found.`

## Warnings

- [ ] Remove or justify `set.seed()` usage in package code.
  - `R/LFQData.R:113`
  - `R/simulate_LFQ_data.R:166`
  - `R/simulate_LFQ_data.R:214`
  - `R/simulate_LFQ_data.R:277`

- [ ] Review deprecated API usage and replace where possible.
  - `.Deprecated()` in `R/LFQData.R:227`
  - `.Deprecated()` in `R/LFQData.R:301`
  - `.Deprecated()` in `R/tidyMS_aggregation.R:826`
  - `.Deprecated()` in `R/tidyMS_aggregation.R:829`

- [ ] Add missing `@return` / value documentation and regenerate Rd files.
  - BiocCheck reports empty or missing `\value` sections across many man pages.
  - Do not edit `man/` directly; update roxygen comments in `R/` and run `make document`.

## Package Metadata Notes

- [ ] Remove `LazyData: true` from `DESCRIPTION` or set it to false.
- [ ] Decide whether to update `Depends: R (>= 4.1)` to Bioconductor's recommended `R (>= 4.5.0)`.
- [ ] Add an `Authors@R` `fnd` role if grant or institutional funding should be credited.
- [ ] Add a `NEWS` file for Bioconductor release announcements.

## Coding Practice Notes

- [ ] Replace `sapply()` with `vapply()` where practical.
  - `R/AnalysisConfiguration.R:266`
  - `R/LFQDataImp.R:95`
  - `R/simulate_LFQ_data.R:141`
  - `R/tidyMS_aggregation.R:712`

- [ ] Replace `1:...` sequences with `seq_len()` or `seq_along()`.
  - `R/LFQDataImp.R:46`
  - `R/squeezeVarRob.R:175`
  - `R/squeezeVarRob.R:452`
  - `R/squeezeVarRob.R:454`
  - `R/tidyMS_missingness_summary.R:253`
  - `R/tidyMS_missingness_summary.R:254`

- [ ] Review `cat()` and `print()` calls outside `show` methods.
  - `R/LFQDataAggregator.R:43`
  - `R/LFQDataPlotter.R:274`
  - `R/LFQDataPlotter.R:307`
  - `R/Model.R:189`
  - `R/ModelFirth.R:194`
  - `R/squeezeVarRob.R:338`
  - `R/squeezeVarRob.R:394`

- [ ] Replace `=` assignment with `<-` where BiocCheck flagged it.
- [ ] Replace `paste()` in condition signals with clearer condition message construction.
- [ ] Remove redundant `stop` / warning calls in condition signals.
- [ ] Avoid direct S4 slot access with `@` in examples and vignettes; use accessors.
- [ ] Remove or justify `<<-` in `R/ContrastsModeratedDEqMS.R:148`.
- [ ] Remove or justify `suppressWarnings()` / message suppression.
  - `R/Contrasts.R:157`
  - `R/tidyMS_build_model.R:267`

## Documentation and Formatting Notes

- [ ] Add runnable examples to exported-object documentation where feasible.
- [ ] Replace `\dontrun{}` with `\donttest{}` in examples.
  - `AggregateLimpa.Rd`
  - `ContrastsLimpaFacade.Rd`
  - Update roxygen sources, not generated `man/` files.
- [ ] Review long functions; BiocCheck reports 38 functions over 50 lines.
  - Longest: `fitFDistRobustly_LG()` in `R/squeezeVarRob.R` at 267 lines.
- [ ] Review Bioconductor formatting notes separately from the package's local 120-character `.lintr` convention.
  - BiocCheck reports 992 lines over 80 characters.
  - BiocCheck reports 60 lines containing tabs.
  - BiocCheck reports 5434 lines whose indentation is not a multiple of 4 spaces.
