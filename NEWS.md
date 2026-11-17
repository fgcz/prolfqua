# prolfqua 1.9.4

## New features

- Fitted analysis results support `print()` for a concise model overview without computing contrasts, and `summary()` for fit status and per-contrast significance counts.

- Fitted analyses retain an independent input snapshot. Matrix models retain their response matrix, weights and observed-response mask, and LOD-refitted models identify imputed rows for result inspection.

- New method `ropecaV2_nested` (`ContrastsROPECAV2NestedFacade`, `ContrastsROPECAV2`): ROPECA with calibrated protein p-values. `ropeca_nested` tests the median signed peptide p-value against a beta distribution, which is conservative for independent peptides and makes proteins with many correlated peptides too significant; its 10-peptide cap only damps that. `ropecaV2_nested` converts the moderated peptide t-statistics to normal scores and tests their median per protein against its exact null distribution, allowing for the correlation of the peptides of a protein, which it estimates once from the peptide model residuals. Share of p < 0.05 for proteins that do not change, by number of peptides (prolfquabenchmark vignette `Benchmark_ropeca`): on simulated proteins whose peptide correlation varies between proteins, 0.050 to 0.067 for 2 to 20 peptides, where `ropeca_nested` gives 0.037 to 0.082; on the human proteins of IonStar, whose peptides are nearly uncorrelated, 0.041 to 0.053 for 2 to more than 20 peptides, where `ropeca_nested` gives 0.018 down to 0.000. `diff` is the median peptide fold change, as in `ropeca_nested`; `statistic`, `std.error` and the 95% interval follow from the p-value. `ropeca_nested` is unchanged.
- New option `prolfqua.verbose`: set `options(prolfqua.verbose = FALSE)` to silence prolfqua's informative messages (columns added, rows removed, files written). The messages carry class `prolfqua_message`, so applications can route them to their own log. Both options are described in `?prolfqua`.

## Behaviour changes

- Source annotations and a standards report now trace all General and Regression rOpenSci standards to evidence, explicit partial limits, non-applicability reasons or deferred work. The feasibility assessment is kept separate from current compliance.

- The package help now defines the API lifecycle, analysis terminology and statistical terms. The Under the hood article explains categorical formula coding and contrast weights, and the analysis help documents measured runtime scaling.

- Invalid analysis, simulation, response, contrast and count inputs now fail at the entry point with a classed error. Method keys must match exactly (case sensitive). A misspelled limma weights column is an error rather than an unweighted fit. Subset scaling rejects mismatched transformation states.

- Fewer hard dependencies: 'UpSetR', 'ComplexHeatmap', 'circlize', 'plotly', 'ggbeeswarm', 'progress', 'logistf', 'lme4' and 'lmerTest' moved from Imports to Suggests, and 'forcats', 'tidyselect' and 'withr' are no longer needed (Imports: 23 to 11). Heatmaps, UpSet plots, interactive plots, the `firth` and `lmer_nested` methods ask to install the packages they need. Without 'plotly' the intensity density plot is static, without 'ggbeeswarm' box plots show jittered points, and without 'progress' no progress bar is shown. Progress of aggregation and of writing box and aggregation plots now also follows `options(prolfqua.progress)`.
- prolfqua no longer prints internal progress messages such as "setup done", "completing cases", "contrasts_linfct", "starting aggregation", the model formula or the number of rows with missing values.
- All errors are raised with `rlang::abort()` and carry class `prolfqua_error`; invalid arguments also carry `prolfqua_error_bad_argument`, so callers can catch prolfqua errors by class. Missing optional packages (`limpa`, `Rfit`) are reported with `rlang::check_installed()`, which offers to install them in interactive sessions.
- `lm_impute`, `rfit_impute`, `limma_impute`, `limma_voom_impute`: the limit-of-detection refit no longer warns "essentially perfect fit" or "non-finite value supplied by optim", and the `limma` variants no longer warn "Partial NA coefficients" for the proteins they then refit. Results are unchanged. A failed per-protein fit is reported as a warning of class `prolfqua_fit_failed`.
- `sim_lfq_data_2factor_config()` no longer warns "Unknown or uninitialised column: `nrPeptides`".
- README: the getting-started example called the removed aggregator method `$medpolish()` and the undeclared package `ggpubr`; it now runs as shown and ends with `build_contrast_analysis()`. Installation is from Bioconductor first.

## Bug fixes

- ANOVA extractors again accept the supported `LRT`, `Rao` and `Cp` tests, with correct p-value column metadata. Direct construction validates the same choices as `get_anova_df()`. Analysis, aggregation and ANOVA selectors reject multi-value arguments instead of silently choosing the first value.
- Statistical standards reports now count only fully met standards as compliant. Partial requirements remain TODOs; fitted-response retention is partial and prediction for new categorical groups is unmet. Internal helper documentation now distinguishes response values, column names, palette sizes and quadrature sizes correctly.
- Condition tests cover all-zero contrasts; covariance-borrowing documentation now describes selection of the largest matching donor subset, and its unreachable error and rescue-skipping branch have been removed.

- The core R6 dependency is now declared explicitly, so subclass-based tests pass dependency checks in fresh environments.

- `deqms` / `deqms_voom` (`moderated_p_deqms()`): proteins left out of the variance-on-count curve, e.g. with zero residual variance, got the prior variance of the first fitted protein instead of the curve's value at their own peptide count, with the warning "'newdata' had 1 row but variables found have N rows". They now get the curve's value. Proteins without residual degrees of freedom get NA moderated statistics without "NaNs produced" warnings; the same holds for the unmoderated statistics of `limma`-based methods.
- `ropeca_nested` / `ContrastsROPECA`: the protein p-value counted peptides without a p-value (failed or non-estimable peptide fits) in the number of peptides, which made those proteins too significant, e.g. p = 0.00067 instead of 0.0031 for two tested peptides and one untested. `summary_ROPECA_median_p_scaled()` now counts only peptides with a p-value in `n_not_na`, as the original ROPECA (PECA) does; `df` of `ropeca_nested` changes with it.

# prolfqua 1.9.3

## Behaviour changes

- The Firth backend fits one model per protein: `build_model_logistf()` requires aggregated data, one row per protein and sample, and stops on peptide-level input. `ModelFirth$models` is the `model_analyse()` result (`model_df`, `strategy`) instead of a list of `models1`, `models2` and `hkey`, and `ContrastsFirth$get_linfct()` returns that one container. `sim_build_models_logistf()` loses its `peptide` argument. `firth` results do not change.
- `ContrastsFirth$get_contrasts()` no longer prints the debug messages "determine linear functions:" and "get_contrasts -> contrasts_linfct".

## Removed

- The `firth_nested` method (`ContrastsFirthNestedFacade`, `build_model_glm_peptide()`). It fitted one 0/1 row per peptide and sample, so its detection p-values became more significant the more peptides a protein had: on the IonStar benchmark, the share of unchanged human proteins with p < 0.05 rose from 1.3% (2-3 peptides) to 7.7% (more than 20). Use `binomial_firth_nested`, which applies the same Firth bias reduction to the detected-peptide counts per protein and sample and takes its degrees of freedom from the samples. `ContrastsHurdleFacade` accepts `binomial_nested` or `binomial_firth_nested` as the detection component.

# prolfqua 1.9.2

## New features

- New `lfqdata_from_wide()`, the inverse of `tidy_to_wide_config()`: one or more feature x sample matrices, their row annotation and their sample annotation back to an `LFQData`. It replaces the matrix-to-long conversions that downstream packages write by hand for AnnData and SummarizedExperiment input.
- `LFQDataTransformer$intensity_array()` takes `transformed`; `intensity_array(exp, force = TRUE, transformed = FALSE)` marks a back-transformation as untransformed, so callers no longer reset `is_transformed(FALSE)` afterwards.
- New `LFQData$is_aggregatable()`: whether there is a hierarchy level below the modelled one to aggregate from.
- `ContrastsFacadeBase$get_missing()` uses the contrast column named in the facade's `ContrastConfiguration` (`contrast_col`) instead of assuming `contrast`, so facades of other packages with their own column name (e.g. SAINT's `Bait`) can inherit it.

## Behaviour changes

- Breaking: a `ContrastsPlotter` histogram spec takes the bin boundaries as `breaks`, e.g. `list(score = "FDR", breaks = seq(0, 1, by = 0.05))`, passed to `ggplot2::geom_histogram()`. The old `xlim = c(0, 1, 0.05)` was not an axis range but a (from, to, by) triple; a spec without numeric `breaks` now stops with `prolfqua_error_bad_argument`.
- Requires tidyr >= 1.3.0: identifier splitting uses `tidyr::separate_wider_delim()` instead of the superseded `separate()`. Results are unchanged; `separate_hierarchy()`, `get_uniprot_id_from_fasta_header()` and the contrast-side parser no longer warn about identifiers with too few or too many pieces.
- The model strategies (`strategy_lm()` and the other `Strategy*` classes), `Model`, `ModelFirth` and `ContrastsFirth` check their arguments when constructed and stop with `prolfqua_error_bad_argument` or `prolfqua_error_missing_column`: a strategy needs a formula, `Model` a model table with the subject columns and a strategy, `ModelFirth` the elements `models1`, `models2` and `hkey`, `ContrastsFirth` a `ModelFirth` and named contrasts. Before, a malformed argument failed later with an unrelated error.
- `setup_analysis()` reports a missing `isotopeLabel` or `qValue` column with a message instead of a warning; both are optional and get the defaults `"light"` and `0`. A missing `nr_children` column still warns, because the default of 1 changes the fitting weights of aggregated data.

## Bug fixes

- `LFQDataToSummarizedExperiment()` stops with an install hint when SummarizedExperiment or S4Vectors is missing; before, it returned `NULL` silently.
- `Contrasts$get_contrasts()` no longer prints the debug lines "determine linear functions:" and "get_contrasts -> contrasts_linfct".
- `firth` and `firth_nested` work when the protein key is not called `protein_Id`: `build_model_logistf()` passed no subject id to `ModelFirth`, so the contrasts looked for a `protein_Id` column and failed.
- Building a contrast analysis (`build_contrast_analysis()` and the facade classes) stops with `prolfqua_error_bad_argument` when a contrast names a group without samples, e.g. `group_A` after all its runs were filtered out. The message names the contrast, the empty group and the groups with samples; before, the analysis went ahead and failed later or returned no estimates for that contrast.
- `LFQDataPlotter$write_pdf()` no longer closes the graphics devices the user has open (it called `graphics.off()` before and after writing). `write_pdf()`, `write_boxplots()`, `LFQDataAggregator$write_plots()` and `LR_test(path = )` close their PDF device when drawing fails, instead of leaving it open so that later plots go into it. `write_boxplots()` returns the file path invisibly.

# prolfqua 1.9.1

## New features

- `LFQData$remove_contaminants()` returns a copy without contaminant proteins, using `config$pattern_contaminants` and the built-in prefixes of `is_contaminant()`, as `remove_decoys()` does for decoys.
- `LFQData$get_Aggregator(method = "limpa")` returns an `AggregateLimpa`, which was not reachable through `get_Aggregator()` before.
- `tidy_to_wide_config()` takes `row_keys`, the columns that identify a row of the wide table; the default is still the hierarchy keys plus the isotope label. Keys that do not identify one value per sample are an error.
- `ContrastsTable$new()` takes an optional `config`, the `ContrastConfiguration` of the stored results.
- `get_contrasts()` takes `all` in every contrast class: `ContrastsROPECANestedFacade$get_contrasts(all = TRUE)` appends the remaining ROPECA columns to the standard ones; `ContrastsHurdleFacade` accepts it and always returns every column.

## Behaviour changes

- Breaking for packages that register their own facades: the `needs` value `"same"` is renamed to `"aggregated"`, so `lookup_facade()` reports `"aggregated"` or `"nested"`, and `register_facade()` rejects any other value with `prolfqua_error_bad_argument`. Before, any string was accepted and every value other than `"nested"` was treated as aggregated input.
- `LFQData$filter_proteins_by_peptide_count()` takes optional `threshold`, `parent` and `child` arguments; the defaults are `config$min_peptides_protein` and the first two hierarchy keys. It now counts the second hierarchy level within the first whatever the `hierarchy_depth`: peptide-level configurations are filtered (before, the call did nothing), and three-level hierarchies modelled at depth 2 count peptides per protein instead of precursors per peptide. A hierarchy with a single level is an error (`prolfqua_error_bad_argument`) instead of a warning.
- `nr_features_experiment()` writes its count to `nr_features_exp` by default instead of `nr_child_exp`, the column of `nr_children_experiment()`, so the two counts can be joined.
- `LFQData$rename_response()` stops with `prolfqua_error_bad_argument` when the new name is another existing column; before, it printed a message and left the data unchanged. Renaming to the current response name does nothing. `rename_response()` and the `is_transformed(value)` setter return the `LFQData` invisibly, so they can be chained.

## Documentation

- The functions that count hierarchy entries (`hierarchy_counts()`, `summarize_hierarchy()`, `HierarchyCountsSample`, `nr_obs_sample()`, `nr_children_experiment()`, `nr_features_experiment()`) link to each other as the "hierarchy counting" family, and each states whether it counts per sample or per experiment, and observed values only or all entries. `nr_obs_sample()` has a title that says what it computes.

# prolfqua 1.9.0

Changes for the Bioconductor review ([BiocContributions #171](https://github.com/Bioconductor/BiocContributions/issues/171)).

## New features

- New method `binomial_firth_nested` (`ContrastsBinomialFirthNestedFacade`, `strategy_binomial_firth()`): the `binomial_nested` detection-count model with Firth's bias reduction in place of the fixed pseudo-count. Like `binomial_nested`, it sums peptide detections per protein and sample, so its dispersion and degrees of freedom come from the samples; unlike `firth_nested`, its p-values therefore do not grow more significant with the number of peptides. Log odds ratios stay finite for proteins detected in one group only.
- New method `hurdle_nested` (`ContrastsHurdleNestedFacade`): a hurdle model from peptide-level input in one call, `build_contrast_analysis(method = "hurdle_nested")`. It aggregates to protein level with median polish, fits `lm_impute` to the protein intensities and `binomial_firth_nested` to the peptide detection, and combines them with `ContrastsHurdleFacade`. Being in the registry, it is also available to prolfquapp as a model.
- `ContrastsHurdleFacade$new(intensity, detection)` combines a fitted intensity facade (`lm`, `lm_impute`, ...) with a fitted `firth_nested`, `binomial_nested` or `binomial_firth_nested` detection facade, as msqrob2's hurdle model does: `p.value` is Fisher's combination of the two component p-values, and `FDR` its BH adjustment within contrast. Proteins only one component estimates keep that component's p-value; fully detected proteins keep their intensity p-value. The component p-values, the detection log odds ratio and the combined statistic are extra columns. No model is fitted.

## Bug fixes

- `citation("prolfqua")` works again: `inst/CITATION` contained a non-ASCII dash that `utils::readCitationFile()` could not read.

## Documentation

- New introductory vignette `vignette("prolfqua")`: the prolfqua hex logo, links to all tutorials on the package website, an introduction centred on `build_contrast_analysis()` with a figure of it as a factory rolling out all modelling methods, installation with BiocManager, and a quick start using `lm_impute`.
- The README and the vignette `ModellingMethods` show the factory figure too; the README links the eight vignettes in reading order, starting with the introduction.
- `?build_contrast_analysis` describes the function as the factory of the modelling methods and lists all 20 built-in methods by family, including `binomial_firth_nested` and `hurdle_nested`, which were missing.
- The package description names `build_contrast_analysis()` and the modelling methods it offers.
- New vignette `vignette("ModellingMethods")`: runs every `build_contrast_analysis()` method on one small simulated dataset, with one summary line per method; it renders in a few seconds.
- The package ships eight vignettes, numbered in reading order: `prolfqua`, `CreatingConfigurations`, `Comparing2Groups`, `ModellingMethods`, `TwoFactorDesigns`, `ContrastsFromFactors`, `ModelBasedImputation` and `QualityControlAndSampleSizeEstimation`; together they build in about 25 seconds. Three more articles (`ComparingMethods`, `SimulateData`, `LowLevelAPI`) are on the website only, <https://fgcz.github.io/prolfqua/articles/>; all keep their URLs.
- The vignettes and articles use `build_contrast_analysis()` throughout. `Comparing2Groups` and "Simulate Peptide level Data" fit their contrasts with `lm_impute` (proteins) and `lm` (peptides) in one call instead of building models and contrasts by hand.
- The articles "Limma Backend", "DEqMS-style Count-Dependent Variance Moderation", "Modelling dataset with two Factors" and "Using limpa Directly" are removed: `ModellingMethods` and the articles "Comparing modelling methods on one dataset" and "Two-factor designs" cover their methods through `build_contrast_analysis()`.
- The imputation article is rewritten and renamed "Filling the gaps: model-based imputation for heatmaps and PCA": why heatmaps and PCA need a complete matrix, why proteins missing in a whole group matter (below the LOD or switched off), how `impute_from_model()` fills the matrix from the `lm_impute` models and the LOD, heatmaps and PCA before and after imputation, and the caveats. It moved to `articles/ModelBasedImputation.html`; the old `TestingMissingInference.html` link redirects there.
- Articles renamed so that file names say what they cover: `ContrastFacades` is now `ComparingMethods` ("Comparing modelling methods on one dataset"), `ContrastFacade2Factor` is `TwoFactorDesigns` ("Two-factor designs"), `Contrasts_from_factors` is `ContrastsFromFactors`. `CreatingConfigurations` is titled "Importing your data: AnalysisConfiguration". The old URLs redirect to the new ones.
- New article "Under the hood: strategy, model and contrasts": the `lm` method built step by step (`strategy_lm()`, `build_model()`, `Contrasts`, `ContrastsModerated`), giving the same result as `build_contrast_analysis()`, with model inspection, unmoderated contrasts and `register_facade()`.
- `Comparing2Groups` describes the simulated data it actually uses (the text still referred to MaxQuant Ionstar data) and no longer looks up protein IDs that do not exist in it. Both vignettes use BiocStyle.
- The vignettes cite the published article (Wolski et al., J. Proteome Res. 2023, 22(4), 1092-1104, doi:10.1021/acs.jproteome.2c00441) instead of the 2022 bioRxiv preprint. `citation("prolfqua")` lists the journal, volume, issue and pages in their own fields, and the README BibTeX no longer swaps volume and issue.
- `?prolfqua` opens a package overview page.
- Typos fixed across the help pages, vignettes and README.

## Behaviour changes

- prolfqua requires R 4.6 or later (was R 4.5.0), matching the Bioconductor release it targets.
- `script_copy_helper_vec()` gains an `overwrite` argument (default `TRUE`, as before).
- `summary_ROPECA_median_p.scaled()` is renamed to `summary_ROPECA_median_p_scaled()` (Bioconductor function names must not contain `.`). There is no alias.

## Removed

- Datasets `data_spectronautDIA250_A`, `data_skylinePRMSample_A` and `data_skylineSRM_HL_A`; their examples use simulated data now. `data_IonstarProtein_subsetNorm` and `data_SAINTe_output` move to prolfquasaint, `data_benchmarkExample` and `data_test_confusion_matrix_scores` to prolfquabenchmark (`prolfqua_data(<name>, package = <package>)`). The unused ExperimentHub manifest `inst/extdata/metadata.csv` is gone. `inst/script/README.md` documents the four datasets that remain.
- The `lm_missing` method (`ContrastsLMMissingFacade`) and `ContrastsMissing` (group-mean substitution with a pooled-variance t-test). `MissingHelpers` keeps `get_stats()` and `get_lod()`; its imputation and contrast methods and the `weighted` argument are gone. All of it lives on in the new package [prolfquatakeoutmissing](https://github.com/prolfqua/prolfquatakeoutmissing): loading it registers `method = "lm_missing"` again. In prolfqua, use `lm_impute`, which refits the proteins plain `lm` cannot estimate with LOD-imputed values.

# prolfqua 1.8.0

This release removes about 3,600 lines (20%) of dead or duplicated code from `R/`. Everything removed had no callers in prolfqua or in the downstream packages (prolfquapp, prophosqua, prolfquappPTMreaders, prolfquasaint, prolfquabenchmark, ptm-pipeline).

## Bug fixes

- `filter_proteins_by_peptide_count()` and `nr_B_in_A()` are back, internal: the `filtered()` method of the saved `data_ionstar` dataset calls the filter, so its removal broke every caller of that dataset, prolfquapp's QC vignette among them. `nr_B_in_A(pdata, level_a, level_b)` only counts the children of each parent; the filter, `LFQData$filter_proteins_by_peptide_count()` and `nr_features_experiment()` count through it instead of each repeating the count.
- `sigma()` for `MASS::rlm` fits took the square root twice. The `rlm` method's contrast `sigma` and all moderated statistics derived from it (`std.error`, `statistic`, `df`, `p.value`, `FDR`, confidence intervals) change; fold changes do not. On simulated data the moderated df drops from about 8,000 to about 25.
- `moderated_p_limma(variance_floor = )` with `robust = FALSE` returned NA moderated p-values for rows above the floor, because the scalar prior df was padded with NA when the floored rows were set to `Inf`. This is the default path of the `binomial_nested` facade (`binomial_bound = TRUE`). Rows above the floor now keep the estimated prior df.
- `merge_contrasts_results()`: when both inputs share a model name, `modelName` is now `<name>_prefer` / `<name>_add` instead of NA for every row.
- `setup_analysis()` no longer fails when `sample_name` and `file_name` are the same column.

## Behaviour changes

- `options(prolfqua.vectorize)` is removed; contrasts always use the per-row path, which was already the default. Setting the option has no effect.
- `ContrastsFirth`: proteins whose contrast fails are reported with the `contrasts_linfct()` warning listing the dropped protein IDs, instead of a message.
- `ContrastsInterface` provides default `get_Plotter()` and `to_wide()` methods, so subclasses no longer need to implement them.
- `estimate_intensity()` takes `method = "medpolish"` or `"rlm"` instead of a function in `.func`.
- A facade's `config` field is filled on the first `get_config()` call; use `get_config()` rather than reading `$config`.
- `center_to_reference_cfg()` always returns a median-centred copy (`summary` and `copy` removed). `scale_with_subset()` gains `colname` and loses the unused `get_scales`.
- `plot_stat_violin_median()` checks `stat` with `match.arg()`, so the default `stat` now works (it used "CV" only when given explicitly).
- `get_p_values_pbeta()` returns NA for an NA `n.obs` instead of stopping.
- `script_copy_helper_vec()` reports the source path it actually copies from and warns once for a missing file; `find_package_file()` returns `NULL` invisibly when the file is not found.

## Removed

- Exported functions: `contrasts_fisher_exact()`, `contrasts_linfct_firth()`, `effective_contaminant_pattern()`, `effective_decoy_pattern()`, `filter_difference()`, `hierarchy_counts_sample()`, `is_singular_lm()`, `lfq_power_t_test_quantiles()`, `linfct_all_possible_contrasts()`, `linfct_factors_contrasts()`, `medpolish_estimate()`, `medpolish_estimate_df()`, `medpolish_estimate_dfconfig()`, `moderated_p_limma_long()`, `names_to_matrix()`, `nr_B_in_A()`, `plot_hierarchies_add_quantline()`, `plot_hierarchies_line_df()`, `plot_lmer_peptide_predictions()`, `remove_na_rows()`, `rlm_estimate()`, `rlm_estimate_dfconfig()`, `scatter_plotly()`, `summarize_stats_all()`, `tidy_to_wide()`, and the function forms of `filter_proteins_by_peptide_count()` and `remove_small_intensities()` (the `LFQData` methods of the same names are unchanged). `unregister_facade()` is no longer exported.
- Replacements where one exists: sample sizes via `LFQDataStats$power_t_test_quantiles()` or `lfq_power_t_test_quantiles_V2()`; `summarize_stats(x, factor_key = NULL)` for `summarize_stats_all(x)`; the Summariser method `$hierarchy_counts_sample()`.
- Datasets `data_checksummarizerobust` and `data_correlatedPeptideList`.
- R6 methods: `LFQData$filter_difference()`, `$omit_na()`, `$contaminant_proportion()`; `LFQDataPlotter$write()`, `$write_pltly()`; `Model$write_anova_figures()`, `$write_coef_figures()` and the same on `ModelFirth`; `ContrastsPlotter$score_plotly()`, `$histogram_estimate()` (use `$histogram_diff()`); `AnalysisConfiguration$id_required()`.
- Unused fields and arguments:
  - strategies: `report_columns` (argument and field of all `strategy_*()` and `Strategy*` classes), `is_mixed`, `$model_fun(get_formula = )` (use `$formula`), `strategy_logistf(test = )`
  - models and contrasts: `ModelLimma$design`, `$formula` and the matching `$new()` arguments; `Model$anova_df`, `ModelFirth$anova_df`; `AnovaExtractor$col_fdr`; `Contrasts$protein_annot`; `ContrastsMissing$global`, `$present`, `$minsd`
  - plotting: `ContrastsPlotter` `prefix` and `protein_annot` fields and `$new(protein_annot = )`; `LFQDataPlotter` `file_paths_pdf`, `file_paths_html`
  - other: `AnalysisConfiguration$hierarchy_keys(rev = )`; `setup_analysis(cc = , from_factors = )`; `plot_hierarchies_line(separate = )`; `transform_work_intensity(intensity_new_name = )`; `compute_pooled(method = )` and `poolvar(method = )` with the "V2" pooled variance; `sim_make_model_lmer(singular = )`; internal facade fields `ContrastsLMMissingFacade$missing_contrast`, `$merged`, `ContrastsLimpaNestedFacade$.lfqdata_nested`
- `build_contrast_analysis()` no longer honours an undocumented `builder` field on facade registry entries.

## Internal

- Facades, strategies, aggregators and model classes now share base classes (`ContrastsFacadeBase$.setup()`, `StrategyBase`, `AggregatorBase`, plotting methods on `ModelInterface`); the four `build_model_limma*()` builders share one core. Class vectors of strategies and aggregators gain `StrategyBase` / `AggregatorBase`.

# prolfqua 1.7.1

- `MissingHelpers$get_lod()` no longer returns NA for data without groups seen once. It falls back to the partly observed groups with the fewest observations, and to the lowest group mean when every group is complete or empty. `build_model_impute()` therefore refits a protein missing a whole group in nearly complete data, where the NA LOD made the refit fail silently and left the protein unfitted.

- New `impute_from_model()` fills the missing response cells of an `LFQData` with predictions from its per-protein `lm` models (from `build_model()` or `build_model_impute()`), leaving observed values unchanged, and reports per protein how many cells were observed and imputed and by which route (`complete`, `fitted`, `lod_refit` or `none`). Additive designs such as `~ group_ + pair` are supported; a model that is not lm-based is an error.

# prolfqua 1.7.0

- A mixed-model analysis (`lmer_nested`) no longer aborts partway through when a
  single protein has no estimable contrast. Such a protein is now dropped with a
  warning naming it, like a rank-deficient one, instead of failing the whole run
  with `subscript out of bounds` from `lmerTest::contest()`.

- A `ContrastConfiguration` built with `pvalue_col = NA_character_` now stores
  `""` instead. HDF5 has no character NA, so a configuration serialized into an
  AnnData artifact came back with `pvalue_col` set to the literal string `"NA"`
  and `has_pvalue()` reported a p-value column that no backend had produced.
  Backends without a p-value, such as SAINTexpress, are now read back
  correctly.

- The peptide count in simulated data is now called `nrPeptides`, matching the
  name every reader in the ecosystem produces. It was the only place the column
  was called `nr_peptides`, which meant code written against simulated data did
  not work on real data -- and vice versa -- and downstream packages carried
  aliases between the two spellings. The `nr_peptides` *option* (minimum
  peptides per protein) is unchanged; only the data column is renamed.

- New `list_to_ContrastConfiguration()` rebuilds a `ContrastConfiguration` from the plain list produced by
  `R6_extract_values()`, so a serialized column-role mapping (e.g. from `SummarizedExperiment` metadata or an AnnData
  `uns` entry) can be restored with its methods, such as `has_pvalue()`, intact.
- Moderated t/Wald contrast tables now report the posterior contrast standard error in `std.error`, while
  `std.error.unmoderated` and `df.unmoderated` preserve the pre-moderation test inputs for downstream comparisons.
  Moderated confidence intervals now use the contrast scale rather than the residual-variance scale; this corrects
  intervals that were too wide when the raw contrast scale factor was below one and too narrow when it was above one.
  The separate infinite-prior-df fallback is unchanged and remains a distinct follow-up decision.
- `main_effect_contrasts()` now emits averaging contrasts with a top-level `LHS - RHS`, making its output directly valid for the stricter contrast parser.
- Added the `binomial_nested` facade for nested LFQ data. It reuses the Firth backend's completed peptide-detection
  representation, fits protein-by-sample detected/undetected peptide counts with a quasibinomial GLM, and reports
  moderated protein-level log odds-ratio contrasts. A symmetric pseudo-count stabilizes separated fits, and posterior
  dispersion is bounded below by one by default (`binomial_bound = FALSE` disables the bound).
- **Variance moderation now uses `limma::squeezeVar()` directly.** `moderated_p_limma()` — and therefore every
  moderated Wald-test facade (`lm`, `rlm`, `rfit`, `lm_impute`, `lm_missing`, `rfit_impute`, `lmer_nested`) and
  the DEqMS small-group fallback — now calls `limma::squeezeVar()` instead of the vendored `squeezeVarRob()` fork
  (limma is already a dependency). Modern limma estimates the empirical-Bayes prior robustly on the fractional / low
  residual degrees of freedom these estimators produce, without the `min_df` workaround prolfqua inherited from MSqRob.
  This slightly changes moderated p-values, FDR and confidence bounds; on the IonStar spike-in benchmark AUC is
  unchanged, FDR calibration is maintained, and ≥99.7% of significance calls are identical. See the
  `SqueezeVar_comparison` vignette in `prolfquabenchmark`. (`ropeca_nested` is unaffected — it aggregates raw
  peptide-level p-values and never used variance moderation. A contrast with fewer than three features now returns
  unmoderated results, matching limma, where the empirical-Bayes prior is undefined.)
- **Fixed the scale used to moderate `rlm` contrasts.** `StrategyRLM$sigma()` now returns the robust scale
  `MASS::rlm` builds its standard errors from (`fit$s`) instead of `stats::sigma()` (the ordinary-residual scale).
  Because the two differ, `ContrastsModerated`'s `sigma / sqrt(var.post)` rescaling was leaving a spurious
  per-protein factor in moderated `rlm` statistics, p-values, FDR and confidence bounds. Only the `rlm` facade's
  moderated output is affected — `lm` and `rfit` were already coherent, and `rlm`'s unmoderated statistics
  (`diff`, `std.error`) were always correct.
- **Removed the exported `squeezeVarRob()`** and its internal helpers (`fitFDist_LG`, `fitFDistRobustly_LG`); the whole
  vendored `squeezeVarRob.R` file is gone. Use `limma::squeezeVar()` instead. DEqMS moderation now calls
  `limma::trigammaInverse()` (identical to the copy that was removed). This is a breaking change only for code that
  called `prolfqua::squeezeVarRob()` directly — the ecosystem does not. Also dropped the now-unused `statmod` from
  `Imports`.
- User-facing errors at prolfqua's public entry points now carry typed condition classes (all inheriting from
  `prolfqua_error`), so callers and tests can catch failures by class instead of matching message text. Covered:
  `setup_analysis()`, `build_contrast_analysis()`, the `LFQData$set_data()` / `set_config_value()` mutators,
  `LFQData$get_Aggregator()`, the aggregators' "must aggregate a multi-level hierarchy" / "aggregate first" checks,
  `LFQDataTransformer$center_to_reference()`, and `ContrastsSimpleImpute$get_contrasts()`. `LFQData$set_data()` and
  `set_config_value()` additionally now reject data that no longer contains the columns required by the current
  configuration, instead of accepting invalid state silently.
- Abundance-density plots now carry per-sample Plotly highlight keys, allowing interactive reports to fade non-hovered
  sample curves while preserving the existing ggplot output.
- **Breaking — contrast schema.** The `modelName` column of `get_contrasts()` output is now the selected facade key (`lm`, `rlm`, `rfit`, `lm_impute`, `lm_missing`, `limma`, `limma_impute`, `limma_voom`, `limma_voom_impute`, `limpa`, `deqms`, `deqms_voom`, `firth`, `lmer_nested`, `ropeca_nested`, `firth_nested`, `limpa_nested`) instead of the testing-schema label (`WaldTest_moderated`, `*_DEqMS`, `*_imputed`, …). Rescue/imputation state moved to a new `estimate_type` column with values `observed`, `lod_imputed`, or `missing_fallback`. The redundant `facade` column was removed (its role is now played by `modelName`). The moderated-Wald-test wording belongs in methods text, not per-row labels. Downstream code that filtered on `modelName == "WaldTest_moderated"` (or the `_imputed`/`_DEqMS` variants) or read the `facade` column must migrate to the facade key and `estimate_type`.
- `build_contrast_analysis()` and the exported `FACADE_REGISTRY` now derive method dispatch and the available-method list from a single seeded registry (`list_facades()`), removing three hand-maintained copies that could drift. New exported base class `ContrastsFacadeBase` holds the shared facade plumbing; the 19 built-in facades are now thin subclasses.
- `ContrastsPlotter` colours volcano/MA/score plots by `estimate_type` when present (and the colour column was left at its default), using black for observed estimates, green for LOD-imputed estimates, and blue for group-mean fallback estimates.
- Correctness fixes: `AggregateTopN` now validates `func` via `match.arg` (an invalid value errors instead of silently meaning `mean`); `Model$get_anova()` drops degenerate rows by `is.na(factor)` rather than a never-matching `"NULL"` string; `sim_lfq_data()` now honours its `mean_prot` argument; `LFQDataTransformer`/`LFQDataSummariser` deep-clone their input like the other decorators (no more mutating the caller's object); `ContrastsFirth$get_linfct()` returns the linfct-annotated model copies instead of mutating the shared model; `contrasts_fisher_exact()` pre-allocates `nrow(x)` results; `ContrastsMissing` now actually validates its output schema.
- Fixed the vectorized contrast path (`options(prolfqua.vectorize = TRUE)`) to match the default loop path on rank-deficient (aliased) model fits: a contrast whose weights fall on non-estimable coefficients now returns `NA` instead of a spurious fold-change of 0 with a `NaN` p-value. The invalid-row guard now counts nonzero weights rather than summing signed weights, so contrasts with canceling `+1`/`-1` weights on missing coefficients are correctly flagged. The default (non-vectorized) path was never affected.
- The vectorized contrast path now assigns stable `contrast_1`, `contrast_2`, ... labels to unnamed contrasts, matching
  the documented naming behavior instead of leaking an expression into the result row names.
- `get_contrast()` now derives `group_1`/`group_2` from the contrast's left/right side expressions, fixing mislabeled per-group columns for averaging contrasts such as `(group_A + group_B)/2 - group_Ctrl` (previously it used the first two group tokens). Simple `A - B` contrasts are unaffected. A contrast that is not a difference `LHS - RHS` now errors with a clear message instead of silently extracting the first two group tokens. Nested contrasts that reference an earlier contrast by name (e.g. `Interaction = "AvsB_gv_X - AvsB_gv_Z"`) remain supported in the `ContrastsMissing` / `lm_missing` path.
- `setup_analysis()` now stops with an informative error (listing the offending keys) when a hierarchy-key/sample combination has more than one observation, instead of silently returning a different-schema count table that crashed downstream. Pass `debug = TRUE` to recover the old behaviour and return the count table for inspection.
- Removed the unused `impute_with_zcomp()`, `estimate_lod_global()`, and `function_lod_quantile()` exports (and the `zCompositions` dependency). For missing-value imputation use `AggregateLimpa$new(lfqdata, impute_only = TRUE)$aggregate()`.
- Hardened `plot_pca()`: errors early on duplicated sample names, an all-missing matrix, or too few samples instead of returning `NULL` (which broke `pca_plotly()`); joins scores to annotation with an explicit `by`; makes `prcomp(center = TRUE, scale. = FALSE)` explicit.
- Hardened abundance heatmaps for sparse significant-feature subsets: when row or column distances are non-finite because of missing values, `plot_heatmap()` now falls back to the input order instead of returning a `ComplexHeatmap` object that fails during drawing.
- Abundance heatmaps now retain every partially observed feature that can be meaningfully row z-scored by default,
  while automatically removing constant, single-observation, and all-missing rows. Significant-feature reports no
  longer hide bait-only results or require a dataset-specific missing-value threshold, and sparse matrices render
  safely even when clustering distances are unavailable.
- `LFQDataPlotter$heatmap()` now shows only the `top_n` most variable features (default 1000), ranked by the prolfqua per-feature statistic (CV for untransformed data, sd for transformed, via `LFQDataStats`). Row clustering uses `stats::hclust`, which errors above 65536 features, so peptide-list / entrapment searches with tens of thousands of degenerate protein groups no longer crash the QC heatmap. Pass `top_n = NULL` (or `Inf`) to keep every feature.
- `StrategyLogistf` now uses Wald confidence intervals instead of profiling every coefficient, preventing `firth_nested`
  analyses from stalling for days on proteins with hundreds or thousands of peptides.

# prolfqua 1.6.1

- Added a `check-bioc` Makefile target and a Bioconductor Docker check image.
- Updated BiocCheck-related documentation and vignette metadata.
- `ContrastsPlotter$volcano()` and `volcano_plotly()` no longer cap positive FDR/p-value scores at `1e-4` by default.
- Renamed `LFQData$to_wide()` to `data_wide()` and removed deprecated compatibility wrappers.
- Split `nr_obs_experiment()` into `nr_children_experiment()` and `nr_features_experiment()`.
- Made `LFQData` data/config internals private and moved callers to accessor methods.
- Refactored transformer and aggregation helpers around current `LFQData` APIs.
