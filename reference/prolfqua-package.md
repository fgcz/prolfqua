# prolfqua: Proteomics Label Free Quantification

Covers the steps of a mass spectrometry based label-free quantification
(LFQ) analysis: quality control, normalization, aggregation of peptides
to proteins, modelling, hypothesis testing and sample size estimation. A
single function, build_contrast_analysis(), fits a model per protein and
computes contrasts with one of many methods: linear, robust and
rank-based regression, 'limma' and 'limma'-voom, 'DEqMS' variance
moderation, the detection probability model of 'limpa', Firth logistic
and binomial models of detection, mixed models and ROPECA on
peptide-level data, a hurdle model combining intensity and detection,
and variants that impute groups without observations at the limit of
detection. All methods return their results through the same interface
and columns, so they can be compared and combined. Models simple designs
with one explanatory variable as well as experiments with several
factors. Developed while analysing dozens of proteomics datasets at the
Functional Genomics Center Zurich.

## Options

- \`prolfqua.verbose\`:

  \`TRUE\` (default) or \`FALSE\`. With \`FALSE\`, prolfqua prints no
  informative messages, such as columns added, rows removed or files
  written. Warnings and errors are not affected. The messages carry
  class \`prolfqua_message\`, so an application can route them to its
  own log with \`withCallingHandlers(prolfqua_message = ...)\`.

- \`prolfqua.progress\`:

  Progress reporting while models and contrasts are computed per
  protein: \`NULL\` (default) shows a terminal progress bar,
  \`"message"\` a throttled message with elapsed time, and a
  \`function(i, total, label)\` is called with the progress.

## Missing and non-finite values

An \`NA\` intensity is a missing observation. The \`method\` argument of
[`build_contrast_analysis`](https://fgcz.github.io/prolfqua/reference/build_contrast_analysis.md)
decides how missing observations are handled: dropped per fit (e.g.
\`lm\`, \`limma\`), imputed at the limit of detection (the \`\_impute\`
methods), or modelled as detection (\`firth\`, \`binomial_nested\`,
\`hurdle_nested\`);
[`impute_from_model`](https://fgcz.github.io/prolfqua/reference/impute_from_model.md)
fills them from fitted models. The \`log\`, \`log2\` and \`log10\`
transforms turn zero intensities (not detected) into \`NA\`. Fitting and
aggregation stop with a classed error on \`Inf\`, \`-Inf\` or \`NaN\`
intensities (\`prolfqua_error_nonfinite_response\`), on \`NA\` in
hierarchy, factor or sample columns (\`prolfqua_error_missing_key\`),
and on inputs without rows (\`prolfqua_error_bad_argument\`). A subject
whose residual variance is zero (constant or noiseless response) keeps
its estimated difference, but its test statistics are \`NA\`, with a
\`prolfqua_zero_variance\` warning, and it is left out of the variance
moderation prior.

## Provenance

Each `LFQData`, model and analysis records the problems and dropped
items that happen in it; `$provenance()` returns them as a table (step,
event, count, IDs, reason). Recorded are values set to \`NA\`
(`remove_small_intensities()`, zeros in log transforms), removed rows,
features and proteins (`remove_decoys()`, `remove_contaminants()`,
[`filter_proteins_by_peptide_count()`](https://fgcz.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md)),
estimation warnings in aggregation (e.g. \`rlm\` not converging),
failed, warned, singular and LOD-refit fits, missing responses dropped
per fit and, for an analysis, the subjects absent from the contrasts
with the reason. Each step reports its new entries once: failures as one
warning of class \`prolfqua_fit_failed\`, everything else as one
\`prolfqua_message\`, which \`prolfqua.verbose = FALSE\` silences.
Derived data (transformed, subset, rows removed) continues the record of
its input; aggregated data and models start their own.

## Lifecycle

prolfqua is stable. Its exported API follows the Bioconductor
deprecation cycle: a renamed or removed interface calls
[`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) for one
release before calling
[`.Defunct()`](https://rdrr.io/r/base/Defunct.html) in a subsequent
release. Changes affecting analyses are recorded in NEWS.

## Glossary

- LFQData:

  A long-format measurement table together with its
  `AnalysisConfiguration`; each row describes a feature in a sample.

- AnalysisConfiguration:

  The mapping from table columns to hierarchy identifiers, sample
  identifiers, categorical factors and the response.

- Hierarchy and hierarchy depth:

  Ordered measurement levels, such as protein, peptide and precursor.
  Depth selects the levels identifying a modelled subject; deeper levels
  identify its child features.

- subject_id:

  The identifier columns of a modelled subject, usually one protein or a
  protein-peptide pair.

- Factor:

  An explanatory variable defining experimental groups.
  [`setup_analysis()`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md)
  converts configured predictors to character labels, including numeric
  and ordered-factor inputs; numeric covariates and ordered-factor
  contrasts are not preserved through this setup.

- Response:

  The numeric intensity column selected by the configuration, usually
  log2 abundance for modelling. On this scale a difference of one
  corresponds to a twofold abundance ratio.

- Contrast:

  A named linear combination of group means or coefficients, for example
  treatment minus control.

- linfct:

  A matrix with one row per requested linear function and one column per
  model coefficient. Multiplying it by the coefficient vector gives the
  estimates being tested.

- Strategy:

  The formula, fitting algorithm and methods for extracting
  coefficients, residual variance, degrees of freedom and contrasts.

- Model:

  Fitted regressions and their diagnostics, either per subject or in a
  matrix representation used by limma.

- Facade:

  An analysis object composing a strategy, model and contrast
  calculation behind the same result and plotting interface.

- Aggregated and nested input:

  Aggregated input has the same hierarchy keys as the output subject.
  Nested input also has child-feature keys; the analysis reports
  parent-level contrasts.

- nr_children:

  The sample-wise number of measured child features, such as peptides
  used to estimate a protein intensity. It may be used as a fitting
  weight; DEqMS separately pools counts per protein.

- LOD and imputation:

  The limit of detection on the working response scale, supplied
  explicitly or estimated from low-observation group means. Imputation
  replaces missing responses. The LOD methods refit otherwise
  unestimable subjects and borrow variance from successful fits.

- estimate_type:

  A result label identifying how an estimate was obtained, for example
  from observed data or an imputed refit.

- Moderation:

  Empirical Bayes shrinkage of residual variances towards a shared
  prior, optionally dependent on abundance or peptide count. Moderated
  degrees of freedom combine residual and prior information.

- Null hypothesis and p-value:

  For a usual two-sided contrast test, the null states that the
  population contrast equals zero. Its p-value is the probability,
  assuming the null and model assumptions hold, of a statistic at least
  as extreme as the observed one. It is not the probability that the
  null hypothesis is true.

- FDR:

  The false discovery rate: the expected proportion of false discoveries
  among rejected hypotheses. The default Benjamini-Hochberg adjustment
  is applied within contrast, assuming independent or suitably
  positively dependent tests. An adjusted p-value is not the probability
  that an individual protein is a false discovery.

- Confidence interval:

  A range from a procedure with the stated repeated-sampling coverage
  under the model assumptions. Wald intervals use the estimated standard
  error and a t quantile; they do not assign a probability to a fixed
  parameter after observing the data.

- Linear-model assumptions:

  On the working scale, classical linear models assume a linear mean,
  independent errors and constant residual variance; exact t inference
  also assumes normal errors. Weights allow specified relative variances
  and mixed models allow specified dependence within subjects. Outliers
  and heavy-tailed errors pull least-squares estimates; the robust
  methods rlm and rfit limit their influence. Robust, detection and
  imputation methods have different assumptions; diagnostics and the
  experimental design must guide the method choice.

## See also

Useful links:

- <https://github.com/fgcz/prolfqua>

- <https://fgcz.github.io/prolfqua/>

- Report bugs at <https://github.com/fgcz/prolfqua/issues>

## Author

**Maintainer**: Christian Panse <cp@fgcz.ethz.ch>
([ORCID](https://orcid.org/0000-0003-1975-3064))

Authors:

- Witold Wolski <witold.wolski@fgcz.uzh.ch>
  ([ORCID](https://orcid.org/0000-0002-6468-120X))

Other contributors:

- Functional Genomics Center Zurich ([ROR](https://ror.org/04tgr1802))
  \[funder\]
