# API Surface Audit: `@export` + `@keywords internal` Functions

**Date:** 2026-03-19
**Scope:** All 113 functions in prolfqua marked with both `@export` and `@keywords internal`
**Downstream packages checked:** prolfquapp, prolfquabenchmark, prolfquappPTMreaders, prolfquasaint, prophosqua

---

## Tier 1: Direct public API — called as `prolfqua::fn()` from downstream packages (11 functions)

These functions are called directly by downstream packages as free functions (`prolfqua::fn()`). They are de facto public API.

**Action taken (2026-03-19):** Removed `@keywords internal` from 9 of 11 functions. The 2 functions called only from `inst/MyArticle/` (`linfct_from_model`, `linfct_matrix_contrasts`) remain internal.

### `R6_extract_values` — AnalysisConfiguration.R

| Package | File | Call |
|---------|------|------|
| prolfquapp | R/R6_AppConfiguration.R:266 | `prolfqua::R6_extract_values(self)` |
| prolfquapp | R/R6_AppConfiguration.R:521 | `prolfqua::R6_extract_values(config)` |
| prolfquapp | inst/application/CMD_MAKE_YAML.R:88 | `prolfqua::R6_extract_values(GRP2)` |
| prolfquapp | inst/application/CMD_DEA_V2.R:145,279,331 | `prolfqua::R6_extract_values(...)` |

### `setup_analysis` — AnalysisConfiguration.R

| Package | File | Call |
|---------|------|------|
| prolfquapp | R/preprocess_BGS_default.R:173 | `prolfqua::setup_analysis(report2, config)` |
| prolfquapp | R/preprocess_DIANN.R:203 | `prolfqua::setup_analysis(peptide, config)` |
| prolfquapp | R/preprocess_MaxQuant.R:421 | `prolfqua::setup_analysis(apeptide, config)` |
| prolfquapp | R/preprocess_FP_PSM.R:577 | `prolfqua::setup_analysis(psma, config)` |
| prolfquapp | R/preprocess_MzMine.R:241 | `prolfqua::setup_analysis(feature, config)` |
| prolfquapp | R/preprocess_CompoundDisc.R:135 | `prolfqua::setup_analysis(peptide, config)` |
| prolfquapp | R/preprocess_MSstats.R:153,262 | Multiple calls |
| prolfquapp | R/proprocess_BGS_default.R:155 | `prolfqua::setup_analysis(report2, config)` |
| prolfquapp | R/tidyMS_MaxQuant.R:439 | `prolfqua::setup_analysis(apeptide, config)` |
| prolfquappPTMreaders | R/preprocess_BGS_site.R:161 | `prolfqua::setup_analysis(site_long, config)` |
| prolfquappPTMreaders | R/preprocess_FP_multisite.R:166 | `prolfqua::setup_analysis(multiSite_long, config)` |
| prolfquappPTMreaders | R/preprocess_FP_combined_STY.R:162 | `prolfqua::setup_analysis(multiSite_long, config)` |
| prolfquasaint | inst/application/SE2_DIANN/DIANN_SE.R:71 | `prolfqua::setup_analysis(peptide, config)` |
| prolfquasaint | inst/application/SE2/CreateSaintExpress_Report.R:88 | `prolfqua::setup_analysis(pdata, config)` |
| prophosqua | vignettes/QCReport.qmd:214 | `prolfqua::setup_analysis(psm, config)` |
| prolfquabenchmark | inst/triqler/Benchmark_Triqler.Rmd:75 | `prolfqua::setup_analysis(resAgr, config)` |
| prolfquabenchmark | vignettes/Benchmark_msqrob2.Rmd, Benchmark_proDA.Rmd, Benchmark_prolfqua.Rmd, etc. | Multiple vignette calls |

### `UpSet_interaction_missing_stats` — tidyMS_missigness.R

| Package | File | Call |
|---------|------|------|
| prolfquapp | inst/application/_Grp2Analysis_V2_metabo_tabs.Rmd:251 | `prolfqua::UpSet_interaction_missing_stats(grp2$RES$lfqData$data, ...)` |
| prolfquapp | inst/application/_Grp2Analysis_V2_R6.Rmd:191 | `prolfqua::UpSet_interaction_missing_stats(dea$lfq_data$data, ...)` |
| prolfquapp | inst/application/_Grp2Analysis_V2_metabo.Rmd:199 | `prolfqua::UpSet_interaction_missing_stats(grp2$RES$lfqData$data, ...)` |
| prolfquapp | inst/application/_Grp2Analysis_Metabolomics.Rmd:181 | `prolfqua::UpSet_interaction_missing_stats(grp2$RES$lfqData$data, ...)` |
| prolfquapp | R/R6_DEAReportGenerator.R:171 | `prolfqua::UpSet_interaction_missing_stats(rd$data, ...)` |
| prolfquasaint | inst/application/SE2_DIANN/DIANN_SE.R:152 | `prolfqua::UpSet_interaction_missing_stats(lfqdataProt$data, ...)` |
| prolfquasaint | inst/application/SE2/CreateSaintExpress_Report.R:140 | `prolfqua::UpSet_interaction_missing_stats(lfqdata$data, lfqdata$config, tr = 2)` |

### `get_UniprotID_from_fasta_header` — utilities.R

| Package | File | Call |
|---------|------|------|
| prolfquapp | R/R6_ProteinAnnotation.R:485 | `prolfqua::get_UniprotID_from_fasta_header(prot_annot, ...)` |
| prolfquasaint | inst/application/SE2/CreateSaintExpress_Report.R:156,162,177 | `prolfqua::get_UniprotID_from_fasta_header(...)` (3 calls) |

### `table_facade` — utilities.R

| Package | File | Call |
|---------|------|------|
| prolfquapp | vignettes/QCandSSE.Rmd:127,242,367,418,435,444 | `prolfqua::table_facade(...)` (6 calls) |
| prolfquabenchmark | inst/triqler/Benchmark_Triqler.Rmd:197,220 | `prolfqua::table_facade(...)` |
| prolfquabenchmark | vignettes/Benchmark_msqrob2.Rmd:245,349,478 | `prolfqua::table_facade(...)` |

### `adjust_p_values` — tidyMS_R6_Modelling.R

| Package | File | Call |
|---------|------|------|
| prolfquabenchmark | vignettes/Benchmark_msqrob2.Rmd:129,443 | `prolfqua::adjust_p_values(all, column = "pval", group_by_col = "contrast")` |

### `matrix_to_tibble` — utilities.R

| Package | File | Call |
|---------|------|------|
| prolfquabenchmark | R/utilities.R:7,8 | `prolfqua::matrix_to_tibble(data.frame(colData))`, `prolfqua::matrix_to_tibble(assay)` |
| prolfquabenchmark | vignettes/Benchmark_msqrob2.Rmd:113,429 | `prolfqua::matrix_to_tibble(hurdle[[i]], preserve_row_names = "name")` |

### `linfct_from_model` — tidyMS_R6_Modelling.R — **kept internal**

Only called from `inst/MyArticle/`, stays `@keywords internal`.

| Package | File | Call |
|---------|------|------|
| prolfquabenchmark | inst/MyArticle/prolfqua_supplement.Rmd:543 | `prolfqua::linfct_from_model(mm$modelDF$linear_model[[1]], as_list = FALSE)` |

### `linfct_matrix_contrasts` — tidyMS_R6_Modelling.R — **kept internal**

Only called from `inst/MyArticle/`, stays `@keywords internal`.

| Package | File | Call |
|---------|------|------|
| prolfquabenchmark | inst/MyArticle/prolfqua_supplement.Rmd:545 | `prolfqua::linfct_matrix_contrasts(linfct, Contrasts)` |

### `pivot_model_contrasts_2_Wide` — tidyMS_R6_Modelling.R

| Package | File | Call |
|---------|------|------|
| prolfquasaint | R/ContrastSaintExpress.R:116 | `prolfqua::pivot_model_contrasts_2_Wide(contrast_minimal, subject_Id = self$subject_Id, ...)` |

### `find_package_file` — vignetteHelpers.R

| Package | File | Call |
|---------|------|------|
| prolfquapp | R/R6_AppConfiguration.R:419 | `prolfqua::find_package_file("prolfquapp", "application/DIANN/config.yaml")` |
| prolfquapp | R/tidyMS_MaxQuant.R:16 | `prolfqua::find_package_file("prolfquapp", "samples/maxquant_txt/tiny2.zip")` |
| prolfquapp | R/preprocess_MaxQuant.R:16 | `prolfqua::find_package_file(...)` |
| prolfquasaint | R/tidyMS_SaintExpress.R:124,126,137,139 | `prolfqua::find_package_file("prolfquasaint", "SaintExpress/bin/...")` (4 calls) |

---

## Tier 2: Internal building blocks — keep `@export` + `@keywords internal` (100 functions)

These functions are used only within prolfqua itself (by R6 classes, other R/ files, vignettes, or tests). They must stay exported because R6 method bodies call them, but they are not direct public API.

This tier includes functions that downstream packages call **only via R6 methods** — the R6 class is the public API, not the underlying free function:

| Function | File | How downstream packages call it |
|----------|------|---------------------------------|
| `complete_cases` | AnalysisConfiguration.R | `lfq$complete_cases()` — R6 method on LFQData |
| `hierarchy_counts` | AnalysisConfiguration.R | `lfqdata$hierarchy_counts()` — R6 method on LFQDataSummariser |
| `hierarchy_counts_sample` | AnalysisConfiguration.R | `sr$hierarchy_counts_sample()`, `sum$plot_hierarchy_counts_sample()` — R6 methods |
| `remove_small_intensities` | tidyMS_R6_TransitionCorrelations.R | `lfqdata$remove_small_intensities()` — R6 method on LFQData |
| `pairs_smooth` | utilities.R | `plotter$pairs_smooth()` — R6 method on LFQDataPlotter |
| `nr_B_in_A` | tidyMS_R6_TransitionCorrelations.R | Only used internally in prolfqua's own `tidyMS_MQ_workflow.R:20` |
| `summary_ROPECA_median_p.scaled` | tidyMS_R6_Modelling.R | Only used internally in prolfqua |

Plus ~93 other functions with no downstream usage at all.

---

## Tier 3: S3 methods (2 functions)

| Function | File | Notes |
|----------|------|-------|
| `sigma.rlm` | tidyMS_R6_Modelling.R | S3 method for `sigma()` on `rlm` objects. Needed for dispatch when `strategy_rlm()` is used. Currently uses `@export`; could be changed to `@exportS3Method`. |
| `df.residual.rlm` | tidyMS_R6_Modelling.R | S3 method for `df.residual()` on `rlm` objects. Same situation. |

---

## Recommended actions (for future implementation)

1. **Tier 1 — done (2026-03-19):** Removed `@keywords internal` from 9 functions called directly as `prolfqua::fn()` by downstream packages. Added missing `@param` tags to pass R CMD check. Two functions (`linfct_from_model`, `linfct_matrix_contrasts`) kept internal — only called from `inst/MyArticle/`.

2. **Tier 2 (100 functions):** Leave as-is. Includes functions reached only via R6 methods — the R6 class is the public interface, the free function is an implementation detail that must stay exported for technical reasons.

3. **Tier 3 (2 S3 methods):** Change `@export` to `@exportS3Method` for proper roxygen2 S3 method registration (not yet done).
