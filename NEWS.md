# prolfqua 1.9.11

## Breaking changes

- `AnalysisConfiguration$ident_q_value` (default column `qValue`, smaller is better) is replaced by `identification_confidence` (default column `identification_confidence`, in [0, 1], larger is better, e.g. 1 - PEP), without backward compatibility: code setting `ident_q_value` must set `identification_confidence` and supply a confidence. `setup_analysis()` sets a missing confidence to 1 instead of a q-value of 0. `data_ionstar` carries `identification_confidence = 1 - pep`.
- Top-N aggregation no longer writes an `ident_qValue` column; every aggregator carries `identification_confidence` instead.

## New features

- New configuration fields `identification_confidence_feature` and `nr_children_feature`: one value per feature, the maximum over samples. `setup_analysis()` derives them when missing, rejects supplied ones with several values per feature, checks that confidences lie in [0, 1] and child counts are non-negative, and reports how many features have a confidence or child count that varies across samples. `complete_cases()` gives added rows their feature's values.
- Aggregation (medpolish, rlm, topN, limpa) gives each protein in each sample the largest identification confidence of its observed peptides, and computes both feature-level columns.
- Models are weighted by identification confidence. Facades that weighted by `nr_children` (`lm`, `lm_impute`, `limma`, `limma_impute`, `limma_voom`, `limma_voom_impute`, `deqms`) weight by `nr_children` times `identification_confidence`; `rlm`, `lmer_nested`, `ropeca_nested`, `ropecaV2_nested`, `deqms_voom`, `limpa` and `limpa_nested` weight by `identification_confidence` and gain a `weights` argument (`NULL` fits unweighted). With confidence 1 everywhere the results are unchanged.
- `weights` of `strategy_lm()`, `strategy_limma()` and the facades accepts several column names, which are multiplied; `strategy_rlm()`, `strategy_lmer()`, `strategy_limpa()` (as `prior.weights`) and `strategy_logistf()` gain a `weights` argument. `sigma()` of an `rlm` fit treats prior weights as inverse variances.
- `ropeca_nested` and `ropecaV2_nested` weight the peptide fits but combine the peptides of a protein without weights. `ContrastsROPECA` and `ContrastsROPECAV2` take `child_weights` to combine them with weights: weighted medians of fold change, statistic and signed p-value or z-score. ROPECA's beta distribution then uses the effective number of peptides (sum w)^2 / sum w^2, rounded down; ROPECA V2 tests the weighted median against its exact null distribution for the given weights (enumerated up to 16 peptides, its normal limit beyond; conservative when a subset of peptides holds exactly half the weight). Under a simulated null with random peptide confidences the weighted combination puts 3.8% (ROPECA) and 7.3% (V2) of p-values below 0.05 and 0.8% and 1.8% below 0.01 (unweighted: 2.5% and 6.3%, 0.6% and 1.6%). In prolfquabenchmark, combining the IonStar MaxQuant peptides weighted by 1 - PEP (0.88 to 1) lowered pAUC_10 by 5.5 (ROPECA) and 4.4 (V2) points; nearly equal weights drop a peptide from ROPECA's effective count and stop V2 averaging the two middle peptides.
- Detection models (`binomial_nested`, `binomial_firth_nested`, and so `hurdle_nested`) count expected true detections: a detected peptide counts with its confidence in the sample, and every peptide counts towards the available ones with its feature-level confidence. `firth` turns each observed protein into a detection with weight c and a non-detection with weight c_feature - c.
- DEqMS moderation takes the peptide count from `nr_children_feature`.
- `rfit` and `rfit_impute` say in a message that they ignore the identification confidence when it is not 1 everywhere; Rfit takes no observation weights.

## Documentation

- The article "Comparing modelling methods on one dataset" fits every method of `build_contrast_analysis()` (all 20, taken from `list_facades()`) and compares them with `ContrastsComparison`: the status of each protein and contrast pair per method, the pairs each method estimates (`$coverage()`, `$upset()`), the pairs of one method (`$pairs()`), the log2 fold changes and statistics on the pairs all intensity methods share, the significant calls, and the sign, significance and coverage conflicts (`$conflicts()`). A last section compares plain limma, run as its User's Guide recommends, as a `ContrastsTable` with its own column names: prolfqua's `limma` without weights reproduces it.

# prolfqua 1.9.10

## Removed

- `ropeca_impute_nested` and `ropecaV2_impute_nested`, added in 1.9.9, are removed without deprecation, because no package used them. In prolfquabenchmark they ranked last or next to last on all three peptide datasets (pAUC_5 by scaled p-value, CPTAC: 25.7 and 37.5, median of all methods 73.7): the peptides refit at the limit of detection weaken the median peptide statistic of truly changed proteins. Use `ropeca_nested` or `ropecaV2_nested`; for proteins without observations in a group, `lm_impute`. The fix of `impute_refit_singular()` from 1.9.9 stays.

# prolfqua 1.9.9

## New features

- New methods `ropeca_impute_nested` and `ropecaV2_impute_nested`: ROPECA and ROPECA V2 on peptide models fitted as in `lm_impute`. A peptide with a group without observations is refit at the limit of detection with a borrowed variance instead of being left out, so more peptides contribute and proteins whose peptides all miss a group get a result; such results have `estimate_type = "lod_imputed"`. V2 estimates the peptide correlation from the observed peptide models only. `ropeca_nested` and `ropecaV2_nested` are unchanged. The factory figure and `facade_styles.tsv` list the 22 methods.

## Bug fixes

- `impute_refit_singular()`, and so `build_model_impute()`, always returns the `imputed` column in the model table; before, it was missing when no fit needed a refit at the limit of detection.

# prolfqua 1.9.8

## New features

- New `ContrastsComparison`: aligns the contrasts of several methods on protein and contrast, from facades or from `ContrastsTable`s holding other tools' results (their `ContrastConfiguration` names the columns). `$coverage()` and `$upset()` show which pairs each method estimates; `$conflicts()` flags pairs where the methods contradict each other at an FDR threshold: significant in opposite directions, significant in one method but not in another that tested the pair, or significant in one and not estimated by another; `$pairs()` returns the pairs of one method with given columns defined. See the section "Comparing methods" of the vignette `ModellingMethods`.

# prolfqua 1.9.7

## Bug fixes

- `rlm` no longer fails with "subscript out of bounds" when its robust weights leave a coefficient of a protein aliased (a regression in 1.9.6, found with the IonStar data of prolfquabenchmark). Contrasts of `lm`, `rlm` and the binomial methods are estimable when they lie in the row space of the weighted design the fit used; contrasts against observations the robust fit weights near zero are not estimable.

# prolfqua 1.9.6

## Behaviour changes

- Proteins with an empty group get exactly the contrasts they can estimate, equal to the differences of their observed group means, in every method. Before, `limma` reported a wrong contrast when the empty group was the reference level (or cell) of the design, e.g. B − Ctrl under the name B − A, with its p-value, and lost the estimable one; `limma_voom` and `deqms_voom` lost every contrast of a protein with any empty group. Subgroup models (`~ subgroup`) were affected for proteins absent from the alphabetically first subgroup. Contrasts a protein cannot estimate have no row in `limma`, `limma_voom`, `deqms_voom`, `lm` and `deqms` (`deqms_voom` and, in `a * b` models, `lm` and `deqms` reported NA rows) and are listed by `get_missing()`. The limma methods (`limma`, `limma_voom`, `deqms_voom`, `limpa` and the impute variants) now compute each protein's contrast standard errors from its own observed, weighted samples; `limma::contrasts.fit()` approximated them for proteins with missing values or sample-varying weights (on simulated data 0.8 to 1.4 times the exact value), so their p-values change. `lm` and `deqms` in `a * b` models with an empty cell now report the cell contrasts that do not need the empty cell (e.g. BZ − BX without AX), which they dropped before.

# prolfqua 1.9.5

## New features

- `inst/extdata/facade_styles.tsv` lists every method of `build_contrast_analysis()` with its family in the factory figure of the README (linear models, limma, limpa, LOD imputation, peptide level, detection, hurdle), the family's colour, and a point shape and line type per method within the family, so plots comparing methods (as in prolfquabenchmark) style each method the same way. The factory figure is drawn from it.
- `LFQDataPlotter$heatmap(row_order = )` draws exactly the given features, in the given order and without row clustering, so two heatmaps of the same features can be compared row by row, e.g. before and after imputation (imputation changes the per-feature sd, so the most variable features and their clustering would differ). `row_order` cannot be combined with `n_most_variable`; an unknown feature ID is an error that names it.

## Behaviour changes

- prolfqua is label-free only: the isotope label is gone. `setup_analysis()` no longer adds an `isotopeLabel` column, so `data_long()`, the wide rowdata and annotation, `hierarchy_counts()`, `summarize_stats()`, the missingness tables and every result table come without it, and `data_wide(as.matrix = TRUE)` row names are the united hierarchy keys (`P1`, not `P1~lfq~light`). An `isotopeLabel` column in the input is an unconfigured column and is dropped; input with two labels per feature and sample fails as duplicate observations. `hierarchy_counts()` returns one row, also for empty input (counts of 0).
- Removed without deprecation, because no package used them: the `AnalysisConfiguration$isotope_label` field, `LFQData$isotope_label()`, and the `isotope_label` argument of `hierarchy_counts()` and `summarize_hierarchy()` (pass `hierarchy` and `factors` by name). `list_to_analysis_configuration()` still loads stored configurations that carry the `isotope_label` key and ignores it.
- prolfqua now requires dplyr >= 1.1.1, ggplot2 >= 3.4.0 and purrr >= 1.0.0; earlier versions lacked functions and arguments prolfqua already used (joins pass `relationship =`, added in dplyr 1.1.1).
- `table_facade(kable = FALSE)` returns the data frame instead of `NULL`. `table_facade()` and `plot_sample_correlation()` (also `LFQDataPlotter$sample_correlation()`) name the Suggests package to install (knitr, corrplot) instead of failing with "there is no package called ...", and every help-page example runs without the Suggests packages.
- `register_facade()` refuses the names of built-in methods (`"lm"`, `"limma"`, ...). Before, registering such a name silently replaced the built-in entry, and the built-in facade then rejected valid input with an error naming the wrong class. Built-in facades check their input shape themselves instead of reading it from the registry.
- Error messages state the problem in their first line ("`method` must be one of ...", "Can't find column `x` in `data`.") and name the function you called instead of prolfqua's internal checker. Condition classes are unchanged; about 30 error and warning messages were reworded. A misspelt `method` lists the methods once and names the value given; before, the list appeared twice.
- Joins inside prolfqua check that their lookup keys are unique: a duplicated subject, contrast or sample key in an input (for instance two rows for one protein and contrast passed to `merge_contrasts_results()`) raises an error instead of silently multiplying rows.
- `list_to_analysis_configuration()` stops with a prolfqua error naming the keys that are not configuration fields (also method names such as `hierarchy_keys`); before, it failed with base R's "cannot add bindings to a locked environment". In the older nested format, keys of the `parameter` block that are not fields (also method names such as `set_response`) are skipped, as other unknown keys there always were; before, a method name failed with "cannot change value of locked binding".
- `moderated_p_deqms_long()` (used by the deqms methods) reports problems of its per-contrast variance fits as a `prolfqua_fit_failed` warning, like the other fitting steps; before, the warning had no class. In the limma methods, limma's warning about proteins without estimable coefficients is their "singular fits" (or "refit at the limit of detection") entry in `$provenance()`, found from the coefficients, not from limma's text; any other limma warning while fitting is a `prolfqua_fit_failed` warning with limma's text as the reason, where before it was a message hidden by `prolfqua.verbose = FALSE`. The variance moderation of `binomial_firth_nested` and `hurdle_nested` reports zero dispersions as a message instead of limma's unclassed "Zero sample variances detected" warning (a warning of class `prolfqua_fit_failed` when more than half are zero); results are unchanged.
- The `deqms` and `deqms_voom` methods and `ContrastsModeratedDEqMS` stop with a `prolfqua_error_constant_counts` error naming the count column when every protein has the same count, e.g. protein data without peptide counts, where `nr_children` is 1; the variance-by-count trend cannot be fitted then. Before, `get_contrasts()` failed inside `predict.loess()` with "NA/NaN/Inf in foreign function call", or returned NaN p-values for every protein.
- Plots: volcano, MA and score plots draw smaller, more transparent points, so dense regions stay readable; `ContrastsPlotter$barplot_threshold()` fills bars with the colour-blind safe palette; score-plot reference lines are black and red regardless of the global `palette()`. `plot_pca()` shapes the second factor with 19 distinguishable shapes and, when that factor has more levels, drops the shape mapping with a message; before, factors with more than 25 levels drew invalid shapes or no points.
- Bad input fails where it enters, with a classed error (`prolfqua_error_bad_argument` unless noted), instead of giving wrong numbers or failing later:
  - `robscale()`, `robscale_subset()` and `get_scales()` stop and name the samples without any measured value; before, one such sample made every scaled value `NA`, without a warning.
  - `LFQData$new()` without setup stops on missing key columns and on more than one row per sample and hierarchy key (`prolfqua_error_duplicate_keys`); before, `get_Stats()` counted the duplicates as replicates.
  - `build_contrast_analysis()` and the facades check the contrasts when the analysis is built. A name that is neither a group (a factor key followed by a level) nor an earlier contrast (e.g. the typo `grp_A`), a contrast whose groups cancel (`group_A - group_A`), a constant, and anything other than `+ - * /` and parentheses (e.g. a function call) are errors. Before, a misspelt group was dropped with a warning and the analysis then failed with an unclassed error, and a contrast string could call any R function.
  - Single-valued arguments are checked: those of the simulators (`n_proteins`, `n_per_group`, `seed`, `with_missing`, ...), `which_missing()`, `LFQData$new()` (`prefix`, `setup`), `setup_analysis(debug)`, `null_split_check()` (`n_splits`, `n_features`, `seed`), the sample-size methods (`delta`, `power`, `sig.level`, `probs`) and the strategies (`model_name`, `trend`, `robust`). Counts must be whole numbers.
  - List columns holding the response or a hierarchy, factor or sample column are refused.
  - `lfqdata_from_wide()` stops when `rowdata` repeats a hierarchy key or the value matrices have different row names; before, repeated keys failed later in `LFQData$new()` with a message about samples, and differing row names went unnoticed.
- `setup_analysis()` keeps factor columns as factors, so their level order sets the reference level of the design and the order of the groups; an ordered factor becomes unordered and unused levels are dropped, both with a message. Before, factor columns became character labels in sorted order. Contrasts do not change; the coefficients of factor input are relative to its first level.
- `LFQData$new()` stores its data as a tibble, also for data.frame and data.table input, and stores the response as double, also when it is an integer.
- Facade constructors name their arguments instead of taking `...`: the limma and limpa methods take `trend` and `robust`, the iterative fits (`rlm`, `rfit`, `rfit_impute`, `lmer_nested`, `binomial_nested`, `binomial_firth_nested`, `firth`) take `control`. `build_contrast_analysis()` passes its `...` to them, so an argument the method does not take is an error; before, `firth` and `hurdle_nested` ignored every further argument, including `control`.
- A Suggests package a method needs (lme4 and lmerTest, logistf, Rfit) is checked once, before any protein is fitted; before, a missing package failed inside the fit of the first protein.
- The `firth` method tests contrasts with the normal distribution (Wald z, `df = Inf`), as its coefficient table, `confint()` and logistf do; before, contrasts used t with n - p degrees of freedom. p-values are smaller, most with few samples; the test stays conservative at small sample sizes. `binomial_firth_nested` and `hurdle_nested`, which moderate an estimated dispersion, keep their t tests.
- The fit status of limpa analyses reports `not_converged` as `NA` instead of 0: limpa returns no convergence status.
- Informative messages name their operation: "Transformed response column added", "filter_proteins_by_peptide_count: column added", "write_boxplots: writing", "write_pdf: writing".
- Plots of data with expected missing values leave them out without ggplot2's "Removed ... rows" warning: the intensity density plot (`LFQDataPlotter$intensity_distribution_density()`), the aggregation plots (`AggregatorBase$plot()`, `plot_estimate()`; lines still break at missing values), the plots of `LFQDataStats` (`violin()`, `density()`, `stdv_vs_mean()`; sd and CV are `NA` for features with fewer than two values) and `missingness_histogram()`.
- Log transforms (`LFQDataTransformer$log2()`, `$intensity_array()`) of negative intensities give one warning of class `prolfqua_negative_intensity` and `NaN`; before, an unclassed warning was followed by R's "NaNs produced".
- `LFQDataPlotter$na_heatmap()` (more than 10,000 rows with missing values) and `LFQDataPlotter$pairs_smooth()` (more than `max` samples) show evenly spaced rows or samples instead of a random selection, so the plot is the same on every call.

## Bug fixes

- `annotation_add_contrasts()` handles designs with more contrasts than samples (e.g. 5 x 2 levels with interactions: 40 contrasts for 30 samples) by appending `NA` sample rows that carry only a contrast, as it already padded the contrast columns for designs with fewer contrasts; before, it failed in `bind_cols()` with "Can't recycle". It stops with a `prolfqua_error_bad_argument` error when two combinations of factor levels get the same group name (`a_b` + `c` and `a` + `b_c` are both `a_b_c`); before, their samples were merged into one group.
- `null_split_check()` no longer adds `pseudo_group_` to the configuration of the `LFQData` passed in; before, later calls on that object such as `data_wide()` and `get_Stats()$stats()` failed.
- `build_model_impute()` stops with an error for input that is not aggregated (subject equals the hierarchy) and for strategies other than `strategy_lm()` and `strategy_rfit()`. Before, peptide-level input was refitted with one placeholder row per missing sample, and mixed models failed with "$ operator not defined for this S4 class".
- `missingness_histogram()` facets by the `factors` it is given; before, it summarised by `factors` but faceted by all configured factors, and failed when `factors` was a subset.
- `table_factors_size()` returns an ungrouped tibble when the design has two or more factors.
- `summarize_stats_quantiles()` (and `LFQDataStats$stats_quantiles()`) return `long` as an ungrouped tibble.
- `summary_ropeca_median_p_scaled()` stops with an error when the `statistic` column is missing; before, it filled the summary's `statistic` column with the column name.
- Hierarchy counts per sample (`HierarchyCountsSample`, `LFQDataSummariser$hierarchy_counts_sample()` and its plot) work without an `nr_children` column, counting every feature once, and stop with an error when `nr_children > 1` is requested without that column.
- `LFQDataPlotter$pca_plotly()` plots the principal components given in `PC`; before, it always plotted the first two.
- `plot_heatmap_cor()` colours by its `color` argument, spread from the lowest to the highest correlation (default white to red, as before); before, `color` was ignored.
- `AggregatorBase$write_plots()` (and the aggregators' `write_plots()`) uses its `show.legend` argument; before, the legend was never shown.
- `linfct_matrix_contrasts()` names an unnamed contrast `contrast_<i>`; before, it dropped it with the warning "LHS of `:=` must be a string".
- Aggregation plots from `plot_estimate()` and `AggregatorBase$plot()`, and boxplots from `LFQDataPlotter$boxplots()`, no longer carry the whole data set: `saveRDS()` of one aggregation plot of 300 proteins writes 11 MB instead of 2 GB.
- `null_split_check(n_splits = 2.5)` stops with an error; before, it ran 98 splits.
- The packaged IonStar peptides (`prolfqua_data("data_ionstar")`) carry the `nr_children` column their configuration names, so they aggregate to proteins; before, aggregation stopped with "Column `nr_children` not found".

## Performance

- The duplicate-key check run before modelling and aggregation (`build_model()`, `build_model_limma()`, the facades and the aggregators), and `print()` of an `LFQData`, count distinct keys by hashing instead of comparing rows; on precursor-level tables with millions of rows this takes a fraction of a second instead of seconds and no longer needs several GB of extra memory.
- `write_plots()` of the aggregators and `LFQDataPlotter$write_boxplots()` build and write one protein at a time, so their memory use no longer grows with the number of proteins (300 proteins: peak 0.33 GB instead of 0.84 GB).
- `summarize_stats()`, and with it `get_Stats()`, the heatmap feature selection and the limit of detection of the impute methods, computes all groups in one pass: 0.2 s instead of 10.5 s for 67,000 peptides in 12 samples.
- Median-polish aggregation is about 4.5 times faster (2,000 proteins: 1.6 s instead of 7.4 s); rlm aggregation fits each protein in base R and is about 3.4 times faster than in the earlier development versions (2,000 proteins: 2.2 s instead of 7.5 s).
- Precursor-scale input (4.8 million rows) no longer makes a copy of the long table in the checks every analysis runs before fitting: the full-rank check of the design takes 0.04 s and no extra memory instead of 2.4 s and 1.4 GB, the check of the contrast groups needs 0.1 GB instead of 0.3 GB, and `nr_obs_sample()` 0.45 GB instead of 0.65 GB. `LFQDataTransformer$robscale()` scales in long format: 0.5 s and 0.25 GB instead of 3.9 s and 1.6 GB. Results are unchanged.
- `LFQDataStats$power_fdr_pilot()` pools the variances of all features at once: 0.1 s instead of 67 s for 100,000 peptides. `Model$get_coefficients()` binds the coefficients of all fits once and is about twice as fast (5,000 proteins: 1.6 s instead of 3.6 s).
- Wald contrasts of the methods that fit one model per protein (lm, rlm, rfit, lm_impute, rfit_impute, deqms, ropeca) are computed once per protein instead of once per contrast row, and contrasts naming a group a protein lacks are skipped without raising and catching an error; the lm method on 5,000 proteins takes 4.9 s instead of 9.4 s.

## Deprecated

- Renamed to snake_case: `lfqdata_to_summarized_experiment()` (was `LFQDataToSummarizedExperiment()`), `lr_test()` (`LR_test()`), `lfq_power_t_test_quantiles()` (`lfq_power_t_test_quantiles_V2()`), `r6_extract_values()` (`R6_extract_values()`), `list_to_analysis_configuration()` (`list_to_AnalysisConfiguration()`), `list_to_contrast_configuration()` (`list_to_ContrastConfiguration()`), `summary_ropeca_median_p_scaled()` (`summary_ROPECA_median_p_scaled()`), and the misspelt `missigness_histogram()` and `LFQDataPlotter$missigness_histogram()` to `missingness_histogram()`. The old names still work and warn; see `?"prolfqua-deprecated"`.
- Simulator arguments renamed: `Nprot` to `n_proteins`, `N` to `n_per_group`, `PEPTIDE` to `peptide`, `TWO` to `two_factor` in `sim_lfq_data()`, `sim_lfq_data_peptide_config()`, `sim_lfq_data_protein_config()`, `sim_lfq_data_2factor_config()` and `sim_build_models_lm()`, `_lmer()`, `_logistf()`. The old arguments still work and warn.
- `volcano_plotly()` and `multigroup_volcano()` take their data as `data`; `.data`, which shadowed the tidy-eval pronoun, still works and warns.
- `merge_contrasts_results(model_name = )` is deprecated: it was never used, and the merged model keeps its name `<prefer>_<add>`. Passing it warns.
- The 73 functions that were exported but documented only on internal help pages are deprecated and will be removed from the exports; prolfqua's interface is `build_contrast_analysis()`, `LFQData` with its `get_*()` objects, and the model and contrast classes. Each still works, warns, and names its public route, e.g. `summarize_stats()`: `LFQData$get_Stats()$stats()`; `plot_pca()`: `LFQData$get_Plotter()$pca()`; `estimate_intensity()`: `LFQData$get_Aggregator("medpolish")$aggregate()`; `compute_contrast()`, `linfct_from_model()`: `build_contrast_analysis()` and `Contrasts$get_linfct()`. The full list, grouped by route, is in `?"prolfqua-deprecated"`. This includes `lfq_power_t_test_quantiles()`, `summary_ropeca_median_p_scaled()` and `missingness_histogram()`, the new names of three functions renamed above; their old names now point to the R6 methods as well.
- `script_copy_helper_vec()` and `separate_hierarchy()` moved to prolfquapp; `get_contrast()` is deprecated without replacement (only prolfquatakeoutmissing used it, which has its own copy); for `matrix_to_tibble()` use `tibble::as_tibble(x, rownames = )`.
- `LFQDataPlotter$heatmap(top_n = )` is renamed `n_most_variable`, since it keeps the features with the largest sd (or CV); `top_n` still works and warns.

## Removed

- The example dataset `x5463yzwer453bbb` (an annotation table of 36 samples). The `annotation_add_contrasts()` example and the article "Contrasts from factors" build their annotation table with `expand.grid()`.
- The old test script in `inst/issue71` is no longer installed.

## Documentation

- `?ContrastsFacadeBase` is a public help page: how a facade for a new modelling backend is written and registered with `register_facade()`, with a working example.
- The reference index of the pkgdown site groups the help pages by topic. The article "Under the hood" explains how contrast strings become weights on the model coefficients with `Contrasts$get_linfct()`, and the imputation vignette draws its before and after heatmaps with `LFQDataPlotter$heatmap(row_order = )`.
- Help pages render Markdown, so code and links that showed as literal backticks or brackets now format; every exported function documents its return value; titles are in sentence case and `lfqdata` arguments share one description.
- Every argument of the exported functions is documented, and the `sim_build_models_logistf()` example checks the fitted models instead of a field that does not exist.
- `?prolfqua` has a section "Input types and conversion": which table classes are accepted, and how the response, factor and hierarchy columns are converted. `?AggregateLimpa` documents the defaults of limpa's iterative steps and why they are not exposed; `?ContrastsFirthFacade` documents its Wald z test. The help pages of the simulators, `LFQData$new()`, `null_split_check()` and the sample-size methods state the type and length of each argument.
- The maintainer is Witold Wolski, the copyright holder the University of Zurich; Christian Panse and Jonas Grossmann are listed as contributors. prolfqua follows the rOpenSci Code of Conduct.
- `?prolfqua` has a Resources section: memory per row, run time per protein for the matrix and per-subject backends, retained fitted models, and the cost of `null_split_check()`. `LFQData$get_config()` says it returns the configuration of the object itself, not a copy, and the ignored `all` argument (kept for the `ContrastsInterface` signature) and `global` argument (kept for compatibility with `Contrasts$get_linfct()`) say so.
- `?build_model_limma_voom` derives the vooma precision weights and the external weights (peptide counts per protein and sample), and states when the fit equals `limma::vooma()` followed by `limma::lmFit()`.
- `ropeca_nested` is documented as a derivation of PECA and ROPECA: `?ContrastsROPECANestedFacade` lists what differs (no ROTS step, any design, at most 10 peptides in the beta correction) and compares it with PECA on the IonStar data: the same fold changes, and the same p-values for proteins with an odd number of peptides up to 10. The script is `data-raw/compare_peca.R`.
- Every public help page has a description, and the return values of the contrast generators, `merge_contrasts_results()`, `pivot_model_contrasts_to_wide()`, `volcano_plotly()` and 15 other functions say what they return. The contrast generator examples and the article "Contrasts from factors" use neutral level names (genotype WT, KO).
- The README is built from `README.Rmd`, so its examples run; it has a statement of need, a list of related packages and the same install command as the vignettes. The vignette "Importing your data" points to the prolfquapp readers for DIA-NN, MaxQuant, Spectronaut and FragPipe output, shows `lfqdata_from_wide()`, and its text matches its code. The reading lists include "Transforming and calibrating LFQ data".
- Vignettes no longer hide warnings. The remaining ones are fits that fail for single proteins (for example `lmer_nested` on proteins with one peptide); the vignettes explain them and show the record in `$provenance()`.
- Corrected help texts: the medpolish `eps` is a relative change of the sum of absolute residuals; `?prolfqua_analysis_methods` says where each backend keeps `converged` and `iterations`.
- The run times in `?prolfqua` (Resources) and `?build_contrast_analysis` are measured again on this version, and the Resources section adds aggregation and the QC writers.
- The GitHub install instructions in the README and the introduction vignette install the remotes package, which `BiocManager::install("fgcz/prolfqua")` needs; before, they failed on an R installation without remotes. The README also says how to check the R version, what to do when prolfqua is not available for the installed Bioconductor version, and that the GitHub version is ahead of the Bioconductor release.

# prolfqua 1.9.4

## New features

- `plot_heatmap_cor()`, `plot_sample_correlation()` and the plotter's `$heatmap_cor()` and `$sample_correlation()` take `method`: a `stats::cor()` method (`"pearson"`, `"spearman"`, `"kendall"`) or a function returning a sample correlation matrix.
- `predict()` on fitted analyses returns predictions for any combination of subject and group (default: every subject in every input sample, including samples with missing responses), on the fitting scale. Groups absent from a subject's fit are predicted from the limit-of-detection refit of the `_impute` methods and are `NA` otherwise; ROPECA predicts per peptide, mixed models at population level, detection models as probabilities, and hurdle analyses per component.
- `LFQData` prints a summary (samples, response and whether it is transformed, hierarchy counts, samples per factor level, response range, missing responses, provenance), and `plot()` on an `LFQData` shows the response distribution per sample.
- `digits` controls the significant digits of numbers printed by `print()` of `LFQData`, of fitted analyses and of their `summary()` (default 3), and of `contrast_summary_table(digits = )`.
- Convergence controls: `strategy_rlm()`, `strategy_lmer()`, `strategy_binomial()`, `strategy_binomial_firth()`, `strategy_logistf()`, `strategy_rfit()`, the `medpolish` and `rlm` aggregators (`get_Aggregator(method, control = )`) and `estimate_intensity()` take a `control` list (iteration limit, tolerance), which `build_contrast_analysis(..., control = )` passes on, also for `firth`. Defaults are unchanged and documented with their reasons; unknown names raise `prolfqua_error_bad_argument`.
- Convergence statistics: model tables keep `converged` and `iterations` per subject (rlm, glm and Firth, lmer, logistf; `NA` where the backend keeps no status), the fit status of `print()` and `summary()` counts `not_converged`, and a fit that stops at its iteration limit without a backend warning is recorded in `$provenance()`.
- Provenance: `LFQData`, models and analyses record problems and dropped items, read with `$provenance()`: values set to `NA`, removed rows, features and proteins, aggregation warnings, failed, warned, singular and LOD-refit fits, missing responses dropped per fit, and subjects absent from the contrasts with the reason. Each step reports its new entries in one message (failures: one warning); `print()` and `summary()` of an analysis show a one-line summary.
- `null_split_check()` checks the calibration of normalization and model together: it splits the samples within each group at random into two pseudo-groups, tests them against each other with any `build_contrast_analysis()` method, and returns per split the share of p < 0.05 and p < 0.01, the number of FDR < 0.05 hits and the Kolmogorov-Smirnov distance from uniform, with the sample assignment and the contrasts. No feature differs between pseudo-groups, so a calibrated pipeline gives about 5% of p < 0.05 and no FDR hits.
- `summarize_variance_explained()` fits a least-squares model per feature of a log-transformed `LFQData` (default: the full factorial of its factors) and returns each term's share of the feature's variance (eta²) and the model R². `pooled_variance_explained()` pools the sums of squares over features; on complete features its R² equals the Euclidean PERMANOVA R². Use them to check that a normalization lowers within-group variance more than overall variance.

- `LFQDataTransformer$subtract_sample_offsets()` applies caller-supplied log2 offsets by sample ID and stores the adjusted intensities in a new response column.

- The Transforming and calibrating LFQ data vignette shows how to use the active response, internal reference proteins and offsets from external calibrants.

- `plot()` now defaults to observed responses and fitted group means with confidence intervals for up to six subjects selected by FDR. Select subjects with `subject` and confidence level with `level`; volcano and MA plots remain available through `type`. Combined hurdle analyses retain the volcano default and expose fitted plots through their intensity and detection components.

- `vcov()` on fitted analyses returns coefficient covariances in a subject-keyed table, using borrowed uncertainty for LOD refits and moderated per-subject covariances for matrix models.

- `fdr_sig_level()` converts a target FDR, null-feature proportion and power into a per-test threshold for the existing protein and SD-quantile sample-size methods. `power_fdr_prior_n()` and `LFQDataStats$power_fdr_pilot()` plan balanced two-group sample sizes by averaging power over a pilot variance prior and a specified effect distribution. The QC and Sample Size Estimation vignette demonstrates both approaches.

- Fitted analyses support `print()` (model overview without computing contrasts), `summary()` (fit status and significance counts per contrast), `coef()`, `confint()` (Wald intervals), `formula()`, `nobs()`, `fitted()` and `residuals()` (with sample identifiers and imputation flags, on the fitted scale) and `plot()` (group means with intervals, volcano or MA plot). All registered backends are tested against this result contract; the Modelling methods vignette shows the methods.

- Fitted analyses retain an independent input snapshot. Matrix models retain their response matrix, weights and observed-response mask, and LOD-refitted models identify imputed rows for result inspection.

- New method `ropecaV2_nested` (`ContrastsROPECAV2NestedFacade`, `ContrastsROPECAV2`): ROPECA with calibrated protein p-values. It converts the moderated peptide t-statistics to normal scores and tests their median per protein against its exact null distribution, allowing for the peptide correlation estimated from the peptide model residuals. `ropeca_nested`, unchanged, is conservative for independent peptides and too significant for proteins with many correlated peptides. Share of p < 0.05 for unchanged proteins (prolfquabenchmark vignette `Benchmark_ropeca`): 0.050 to 0.067 on simulated data and 0.041 to 0.053 on IonStar, against 0.037 to 0.082 and 0.018 to 0.000 for `ropeca_nested`. `diff` is the median peptide fold change; `statistic`, `std.error` and the 95% interval follow from the p-value.
- New option `prolfqua.verbose`: set `options(prolfqua.verbose = FALSE)` to silence prolfqua's informative messages (columns added, rows removed, files written). The messages carry class `prolfqua_message`, so applications can route them to their own log. Both options are described in `?prolfqua`.

## Behaviour changes

- `LFQDataStats$power_t_test()`, `$power_t_test_quantiles()` and `$power_fdr_pilot()` raise a `prolfqua_error_bad_argument` error on untransformed data. Before, they warned "...; aborting." and returned `NULL`, which a pipeline could pass on unnoticed.
- All ggplot figures have a new default look in the style of GraphPad Prism: white background, black axis lines with outward ticks, bold axis titles, bold facet labels on a light grey strip. Volcano, MA, score, fitted-value, PCA, density and SD-versus-mean plots keep faint gridlines; histograms, bar charts, boxplots and violins have none. Violins have a thinner outline (`linewidth = 0.3`). Add any ggplot2 theme to override, or add the new `theme_prolfqua()` to your own ggplots to match.
- `setup_analysis()` announces each transformation it applies with a `prolfqua_message`: source columns united into one key, non-character factor columns converted to character labels, unconfigured columns dropped, and rows with NA response added for absent feature and sample combinations. Its help lists every step and how to avoid it.
- `build_contrast_analysis()` checks before fitting that the fixed-effect design is of full rank across samples (aliased factors, empty interaction cells) and stops with `prolfqua_error_collinear_design` naming the aliased coefficients; such designs used to fail later with unrelated errors.
- Fitting and aggregation reject data with more than one row per sample and hierarchy key (`prolfqua_error_duplicate_keys`), also when the `LFQData` was created without `setup_analysis()`.
- `LFQData$new()` drops extra classes and attributes (e.g. `I()`, units, labels) from the response column with a message; transformations and models use the plain numeric values.
- Colour-blind safe default colours: abundance heatmaps use blue-white-red instead of green-black-red, factor levels in heatmap annotations and PCA use Okabe-Ito (viridis beyond 8 levels) instead of the ggplot2 hue palette, estimate types in volcano, MA and score plots use black, orange and blue instead of black, green and blue, and zero lines are grey instead of green.
- Larger default text: `multigroup_volcano()` labels default to `size = 3` (was 1); `plot_sample_correlation()` uses the corrplot text sizes and passes `...` (e.g. `tl.cex`, `number.cex`) to `corrplot::corrplot.mixed()`.
- `?prolfqua` and `?LFQData` state how violated assumptions affect estimates: outliers pull least-squares fits while `rlm` and `rfit` resist them, and missing low intensities shrink fold changes of down-regulated proteins toward zero, also with LOD imputation.

- Analysis help and the Modelling methods vignette document internal aggregation, detection encoding, retained responses and output scales, including which steps have no inverse.

- The package help now defines the API lifecycle, analysis terminology and statistical terms. The Under the hood article explains categorical formula coding and contrast weights, and the analysis help documents measured runtime scaling.

- Invalid analysis, simulation, response, contrast and count inputs now fail at the entry point with a classed error. Method keys must match exactly (case sensitive). A misspelled limma weights column is an error rather than an unweighted fit. Subset scaling rejects mismatched transformation states.

- Fewer hard dependencies (Imports: 23 to 12): 'UpSetR', 'ComplexHeatmap', 'circlize', 'plotly', 'ggbeeswarm', 'progress', 'logistf', 'lme4' and 'lmerTest' moved to Suggests; 'forcats' and 'tidyselect' are no longer needed. Functions that need a suggested package ask to install it; without 'plotly', 'ggbeeswarm' or 'progress' the density plot is static, box plots show jittered points and no progress bar is shown. Progress of aggregation and plot writing follows `options(prolfqua.progress)`.
- prolfqua no longer prints internal progress messages such as "setup done", "completing cases", "contrasts_linfct", "starting aggregation", the model formula or the number of rows with missing values.
- All errors are raised with `rlang::abort()` and carry class `prolfqua_error`; invalid arguments also carry `prolfqua_error_bad_argument`, so callers can catch prolfqua errors by class. Missing optional packages (`limpa`, `Rfit`) are reported with `rlang::check_installed()`, which offers to install them in interactive sessions.
- `lm_impute`, `rfit_impute`, `limma_impute`, `limma_voom_impute`: the limit-of-detection refit no longer warns "essentially perfect fit" or "non-finite value supplied by optim", and the `limma` variants no longer warn "Partial NA coefficients" for the proteins they then refit. Results are unchanged.
- Missing and non-finite values: log transforms turn zero intensities into `NA` with a message instead of `-Inf` with a warning. Fitting and aggregation stop with classed errors on `Inf`, `-Inf` or `NaN` intensities (`prolfqua_error_nonfinite_response`), on `NA` in hierarchy, factor or sample columns (`prolfqua_error_missing_key`) and on inputs without rows, where before the protein vanished from the contrasts behind a warning. `which_missing()` rejects input containing `NA`. See `?prolfqua` and the Modelling methods vignette.
- Zero residual variance: a subject with a constant or noiseless response (e.g. all intensities equal, or identical within each group) keeps its `diff`, but its standard error, statistic, p-value and confidence limits are `NA`, with one `prolfqua_zero_variance` warning per contrast table. Before, these tests divided floating-point noise by noise (an rlm p-value of 1e-17 for a constant protein), and the near-zero variances shrank the moderation prior of all other proteins. Applies to lm, rlm, rfit, lmer and the limma and DEqMS methods; detection models keep their dispersion floor.
- `rlm` on a constant or noiseless subject returns the least-squares estimates with robust scale 0 instead of iterating `MASS::rlm()`, whose scale could underflow to zero and turn the estimates `NA` on Linux and Windows.
- Failed fits and warnings while fitting or aggregating (e.g. non-converging `rlm`, `lmer`, Firth or `medpolish`) are reported once per step as one `prolfqua_fit_failed` warning naming subjects and reasons, instead of one warning per protein; lmer contrast warnings join the single `contrasts_linfct` warning, and limma's "Partial NA coefficients" warning becomes a "singular fits" provenance entry. `transform_work_intensity()` returns the rows whose zeros became `NA` (`zeros`) and `estimate_intensity()` the per-subject warnings (`failures`) instead of signalling them.
- `sim_lfq_data_2factor_config()` no longer warns "Unknown or uninitialised column: `nrPeptides`".
- README: the getting-started example called the removed aggregator method `$medpolish()` and the undeclared package `ggpubr`; it now runs as shown and ends with `build_contrast_analysis()`. Installation is from Bioconductor first.

## Bug fixes

- `adjust_p_values()`, `summarize_hierarchy()`, `missingness_per_condition()`, `missingness_per_condition_cumsum()`, `HierarchyCountsSample` and the `summary` tables of `ContrastsPlotter$barplot_threshold()` return ungrouped tibbles. Before, they kept the grouping used inside, so a later `mutate()` or `summarise()` on the result ran per group without asking for it.
- `LFQData$get_subset()` errors when `x` has none of the subject-id columns. Before, it joined on no columns and returned the whole dataset.
- A data column named `threshold` or `N` no longer replaces the argument of the same name in `remove_small_intensities()`, `filter_proteins_by_peptide_count()` and top-N aggregation; the same holds for `limit` in `LFQDataPlotter$pairs_smooth()`.
- Plots: `ContrastsPlotter$histogram_diff()` marks each contrast's own median (before: the median over all contrasts, drawn once per row); volcano, MA and score plots accept a colour column with more than 4 levels (before: "Insufficient values in manual scale") and `multigroup_volcano()` accepts significance column names that are not syntactic, such as `adj p`; `plot_stat_violin_median()` shows median and quartiles as one point range (before: three unlabelled points per group); beeswarm points in `plot_hierarchies_boxplot()` are 1 to 4 mm (before: the child count in mm).
- `medpolish` and `rlm` aggregation no longer report "Response `medpolish`: dropped attributes names" (or `lmrob`): the estimates are returned as plain numeric values.
- `rlm` aggregation (`get_Aggregator("rlm")`) removed estimable design columns when a sample or peptide column was a factor with an unused level, and then silently used per-sample means for that protein. It now drops unused levels and removes only aliased columns. Data set up with `setup_analysis()` (character columns) is not affected.
- `LFQDataPlotter$write_pdf()` and `$write_boxplots()` no longer append `.pdf` to a file name that already ends in `.pdf`.
- `get_complete_model_fit()`, and with it `Model$get_anova()`, `Model$get_coefficients()` and the coefficient and ANOVA plots, now leaves out fits with fewer coefficients than the full design, for example a protein with one group entirely missing. Before, every fit counted as complete, so these tables also listed such proteins, whose coefficients referred to another reference level. Limit-of-detection imputation is unchanged.

- `moderated_p_limma()`: when all variances are equal, one missing residual df no longer makes the prior df, and with it every moderated p-value, `NA`; the fallback prior df counts only rows with a variance.
- `sim_lfq_data_peptide_config()`, `sim_lfq_data_protein_config()` and `sim_lfq_data_2factor_config()` with `with_missing = FALSE` return complete data with positive peptide counts; before, intensities chosen as missing were kept with a peptide count of 0, which gave them fitting weight 0.

- The introductory and Modelling methods vignettes no longer fail when `R CMD check --no-build-vignettes` runs their code: the factory figure is left out of the extracted R code. The Transforming and calibrating LFQ data vignette ends with session information.

- Coefficient tables for LOD-refitted proteins use their borrowed covariance and degrees of freedom. Robust and rank-based models expose consistent coefficient and p-value columns, restoring their coefficient volcano plots. Firth coefficient tables skip failed fits.

- Printing fitted analyses extracts coefficients from only the first few subject fits for its preview. Limma summaries correctly count failed fits and report singularity as unknown (`NA`) instead of inferring it from failed fits.

- ANOVA extractors again accept the supported `LRT`, `Rao` and `Cp` tests, with correct p-value column metadata. Direct construction validates the same choices as `get_anova_df()`. Analysis, aggregation and ANOVA selectors reject multi-value arguments instead of silently choosing the first value.
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
- Class vectors of strategies and aggregators include `StrategyBase` and `AggregatorBase`.

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
