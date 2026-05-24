# Contrast Configuration

Holds the column-role mapping and behaviour flags that describe how a
specific contrast backend (lm, limma, SAINT, ...) lays out its result
table. Mirrors
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md)
for the modelling side: it carries names, not data.

Consumers that need to operate on a contrast table generically (filter,
build ORA inputs, build rank lists, summarise per protein) should read
column names through the accessors on this object rather than
hard-coding them. This lets backends with different column conventions
(e.g. SAINT's `BFDR`/`log2_EFCs`) plug into the same downstream code
path as LM-style backends without renaming columns.

Defaults match the standard prolfqua contrast schema produced by
[`Contrasts`](https://wolski.github.io/prolfqua/reference/Contrasts.md),
[`ContrastsLimma`](https://wolski.github.io/prolfqua/reference/ContrastsLimma.md),
[`ContrastsModerated`](https://wolski.github.io/prolfqua/reference/ContrastsModerated.md)
etc.

## Value

An R6 class generator.

## See also

Other configuration:
[`AnalysisConfiguration`](https://wolski.github.io/prolfqua/reference/AnalysisConfiguration.md),
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`R6_extract_values()`](https://wolski.github.io/prolfqua/reference/R6_extract_values.md),
[`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md),
[`make_interaction_column()`](https://wolski.github.io/prolfqua/reference/make_interaction_column.md),
[`make_reduced_hierarchy_config()`](https://wolski.github.io/prolfqua/reference/make_reduced_hierarchy_config.md),
[`sample_subset()`](https://wolski.github.io/prolfqua/reference/sample_subset.md),
[`separate_hierarchy()`](https://wolski.github.io/prolfqua/reference/separate_hierarchy.md),
[`setup_analysis()`](https://wolski.github.io/prolfqua/reference/setup_analysis.md),
[`table_factors()`](https://wolski.github.io/prolfqua/reference/table_factors.md),
[`table_factors_size()`](https://wolski.github.io/prolfqua/reference/table_factors_size.md)

## Public fields

- `subject_id`:

  columns identifying the contrast subject (e.g. protein id)

- `model_name_col`:

  column carrying the model name label

- `contrast_col`:

  column carrying the contrast label

- `effect_col`:

  column with the (signed) effect size; for LM this is `diff`, for SAINT
  `log2_EFCs`

- `score_col`:

  column with the per-contrast test statistic

- `pvalue_col`:

  column with the raw p-value, or `NA_character_` for backends that do
  not produce one (e.g. SAINTexpress)

- `fdr_col`:

  column with the FDR / BFDR / adjusted p-value

- `avg_abundance_col`:

  column with the per-feature mean abundance

- `supports_dea_qc`:

  whether the backend supports the legacy differential-expression QC
  HTML report

- `needs_saint_annotation`:

  whether the backend needs annotation reading in SAINT mode (bait /
  control columns)

- `significance_directional`:

  if `TRUE`, the backend only considers positive-effect features
  biologically meaningful (e.g. SAINTexpress, where negative `log2_EFCs`
  are not interpreted as "down" hits). Affects `filter_significant()`:
  when `TRUE` the filter is one-sided on the effect column; when `FALSE`
  it is symmetric (`|effect| > threshold`).

## Methods

### Public methods

- [`ContrastConfiguration$new()`](#method-ContrastConfiguration-new)

- [`ContrastConfiguration$has_pvalue()`](#method-ContrastConfiguration-has_pvalue)

- [`ContrastConfiguration$clone()`](#method-ContrastConfiguration-clone)

------------------------------------------------------------------------

### Method `new()`

Construct a ContrastConfiguration. All arguments map directly to fields
of the same name.

#### Usage

    ContrastConfiguration$new(
      subject_id = character(),
      model_name_col = "modelName",
      contrast_col = "contrast",
      effect_col = "diff",
      score_col = "statistic",
      pvalue_col = "p.value",
      fdr_col = "FDR",
      avg_abundance_col = "avgAbd",
      supports_dea_qc = TRUE,
      needs_saint_annotation = FALSE,
      significance_directional = FALSE
    )

#### Arguments

- `subject_id`:

  columns identifying the contrast subject

- `model_name_col`:

  model-name column

- `contrast_col`:

  contrast-label column

- `effect_col`:

  signed effect column

- `score_col`:

  test-statistic column

- `pvalue_col`:

  raw p-value column, or `NA_character_`

- `fdr_col`:

  FDR / adjusted-p column

- `avg_abundance_col`:

  mean-abundance column

- `supports_dea_qc`:

  supports DEA QC HTML

- `needs_saint_annotation`:

  needs SAINT-mode annotation

- `significance_directional`:

  one-sided filter on effect column

------------------------------------------------------------------------

### Method `has_pvalue()`

Does the backend produce a usable p-value column?

#### Usage

    ContrastConfiguration$has_pvalue()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ContrastConfiguration$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
cfg <- ContrastConfiguration$new(subject_id = "protein_Id")
cfg$effect_col
#> [1] "diff"
cfg$fdr_col
#> [1] "FDR"
cfg$has_pvalue()
#> [1] TRUE
cfg$supports_dea_qc
#> [1] TRUE

# SAINT-flavoured config
saint <- ContrastConfiguration$new(
  subject_id = "protein_Id",
  contrast_col = "Bait",
  effect_col = "log2_EFCs",
  score_col = "SaintScore",
  pvalue_col = NA_character_,
  fdr_col = "BFDR",
  supports_dea_qc = FALSE,
  needs_saint_annotation = TRUE
)
saint$has_pvalue()
#> [1] FALSE
```
