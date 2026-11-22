# LFQData R6 class

LFQData R6 class

## Value

An R6 class generator.

## Missing Data Assumptions

The filtering and imputation methods in this package assume that missing
values are Missing Completely At Random (MCAR) or Missing At Random
(MAR). Abundance-dependent missingness (MNAR), which is common in DDA
proteomics, is not modelled. When low intensities are missing, fold
changes of down-regulated proteins shrink toward zero; LOD imputation
(the `_impute` methods) caps estimates at the LOD and does not remove
this bias. MNAR can also inflate false discovery rates. How `NA`, zero
and non-finite intensities are treated is described in
[`prolfqua-package`](https://fgcz.github.io/prolfqua/reference/prolfqua-package.md),
section "Missing and non-finite values".

## See also

Other LFQData:
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md),
[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md),
[`LFQDataPlotter`](https://fgcz.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md),
[`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md),
[`lfqdata_to_summarized_experiment()`](https://fgcz.github.io/prolfqua/reference/lfqdata_to_summarized_experiment.md),
[`plot.LFQData()`](https://fgcz.github.io/prolfqua/reference/plot.LFQData.md)

## Public fields

- `prefix`:

  e.g. "peptide\_", "protein\_", "compound\_"

## Methods

### Public methods

- [`LFQData$new()`](#method-LFQData-initialize)

- [`LFQData$provenance()`](#method-LFQData-provenance)

- [`LFQData$get_provenance_log()`](#method-LFQData-get_provenance_log)

- [`LFQData$print()`](#method-LFQData-print)

- [`LFQData$set_data()`](#method-LFQData-set_data)

- [`LFQData$get_config()`](#method-LFQData-get_config)

- [`LFQData$set_config_value()`](#method-LFQData-set_config_value)

- [`LFQData$get_copy()`](#method-LFQData-get_copy)

- [`LFQData$is_aggregatable()`](#method-LFQData-is_aggregatable)

- [`LFQData$get_sample()`](#method-LFQData-get_sample)

- [`LFQData$get_subset()`](#method-LFQData-get_subset)

- [`LFQData$subject_id()`](#method-LFQData-subject_id)

- [`LFQData$is_transformed()`](#method-LFQData-is_transformed)

- [`LFQData$remove_small_intensities()`](#method-LFQData-remove_small_intensities)

- [`LFQData$filter_proteins_by_peptide_count()`](#method-LFQData-filter_proteins_by_peptide_count)

- [`LFQData$remove_decoys()`](#method-LFQData-remove_decoys)

- [`LFQData$remove_contaminants()`](#method-LFQData-remove_contaminants)

- [`LFQData$decoy_proportion()`](#method-LFQData-decoy_proportion)

- [`LFQData$complete_cases()`](#method-LFQData-complete_cases)

- [`LFQData$data_wide()`](#method-LFQData-data_wide)

- [`LFQData$factors()`](#method-LFQData-factors)

- [`LFQData$hierarchy()`](#method-LFQData-hierarchy)

- [`LFQData$response()`](#method-LFQData-response)

- [`LFQData$hierarchy_keys()`](#method-LFQData-hierarchy_keys)

- [`LFQData$relevant_hierarchy_keys()`](#method-LFQData-relevant_hierarchy_keys)

- [`LFQData$factor_keys()`](#method-LFQData-factor_keys)

- [`LFQData$relevant_factor_keys()`](#method-LFQData-relevant_factor_keys)

- [`LFQData$sample_name()`](#method-LFQData-sample_name)

- [`LFQData$file_name()`](#method-LFQData-file_name)

- [`LFQData$nr_children_col()`](#method-LFQData-nr_children_col)

- [`LFQData$isotope_label()`](#method-LFQData-isotope_label)

- [`LFQData$data_long()`](#method-LFQData-data_long)

- [`LFQData$rename_response()`](#method-LFQData-rename_response)

- [`LFQData$hierarchy_counts()`](#method-LFQData-hierarchy_counts)

- [`LFQData$summarize_hierarchy()`](#method-LFQData-summarize_hierarchy)

- [`LFQData$get_Plotter()`](#method-LFQData-get_Plotter)

- [`LFQData$get_Summariser()`](#method-LFQData-get_Summariser)

- [`LFQData$get_Stats()`](#method-LFQData-get_Stats)

- [`LFQData$get_Transformer()`](#method-LFQData-get_Transformer)

- [`LFQData$get_Aggregator()`](#method-LFQData-get_Aggregator)

- [`LFQData$clone()`](#method-LFQData-clone)

------------------------------------------------------------------------

### `LFQData$new()`

initialize

#### Usage

    LFQData$new(data, config, prefix = "ms_", setup = FALSE, provenance = NULL)

#### Arguments

- `data`:

  data.frame

- `config`:

  configuration

- `prefix`:

  will be use as output prefix

- `setup`:

  is data setup needed, default = FALSE, if TRUE, calls
  [`setup_analysis`](https://fgcz.github.io/prolfqua/reference/setup_analysis.md)
  on data first.

- `provenance`:

  internal: record of an LFQData this one continues (the same data
  transformed or subset); a copy is kept. NULL starts a new record.

------------------------------------------------------------------------

### `LFQData$provenance()`

problems and dropped items recorded for this data: values set to NA,
rows, features or proteins removed

#### Usage

    LFQData$provenance()

#### Returns

tibble with columns step, event, n, ids (list), reason, failure,
response

------------------------------------------------------------------------

### `LFQData$get_provenance_log()`

internal: the record this data continues, passed to `LFQData$new()` for
derived data

#### Usage

    LFQData$get_provenance_log()

#### Returns

Provenance object

------------------------------------------------------------------------

### `LFQData$print()`

print samples, response, hierarchy counts, samples per factor level,
missing responses and provenance

#### Usage

    LFQData$print(digits = 3, ...)

#### Arguments

- `digits`:

  significant digits of printed numbers (default 3)

- `...`:

  unused

#### Returns

self (invisible)

------------------------------------------------------------------------

### `LFQData$set_data()`

set data (replaces the internal data frame)

#### Usage

    LFQData$set_data(new_data)

#### Arguments

- `new_data`:

  data.frame

#### Returns

self (invisible)

------------------------------------------------------------------------

### `LFQData$get_config()`

return the AnalysisConfiguration object

#### Usage

    LFQData$get_config()

#### Returns

the AnalysisConfiguration of this object itself, not a copy: changing a
field of it changes this `LFQData`. Use
`get_config()$clone(deep = TRUE)` to modify a configuration for another
object.

------------------------------------------------------------------------

### `LFQData$set_config_value()`

set a config field value

#### Usage

    LFQData$set_config_value(field, value)

#### Arguments

- `field`:

  character - field name

- `value`:

  the value to set

------------------------------------------------------------------------

### `LFQData$get_copy()`

get deep copy

#### Usage

    LFQData$get_copy()

------------------------------------------------------------------------

### `LFQData$is_aggregatable()`

can the data be aggregated, i.e. is there a hierarchy level below the
modelled one (`hierarchy_depth`)

#### Usage

    LFQData$is_aggregatable()

#### Returns

logical

------------------------------------------------------------------------

### `LFQData$get_sample()`

samples subset of data

#### Usage

    LFQData$get_sample(size = 100, seed = NULL)

#### Arguments

- `size`:

  size of subset default 100

- `seed`:

  set seed

------------------------------------------------------------------------

### `LFQData$get_subset()`

get subset of data

#### Usage

    LFQData$get_subset(x)

#### Arguments

- `x`:

  data frame with columns containing subject_id

------------------------------------------------------------------------

### `LFQData$subject_id()`

get subject ID columns

#### Usage

    LFQData$subject_id()

------------------------------------------------------------------------

### `LFQData$is_transformed()`

is data transformed

#### Usage

    LFQData$is_transformed(is_transformed)

#### Arguments

- `is_transformed`:

  logical

#### Returns

logical when called without argument, otherwise self (invisible)

------------------------------------------------------------------------

### `LFQData$remove_small_intensities()`

some software is reporting NA's as 0, you must remove it from your data

#### Usage

    LFQData$remove_small_intensities(threshold = 4)

#### Arguments

- `threshold`:

  default 4.

#### Returns

self

------------------------------------------------------------------------

### `LFQData$filter_proteins_by_peptide_count()`

remove parents (e.g. proteins) with fewer than `threshold` distinct
children (e.g. peptides) and add the count as a column. By default the
second hierarchy level is counted within the first, independent of
`hierarchy_depth`; a hierarchy with one level is an error.

#### Usage

    LFQData$filter_proteins_by_peptide_count(
      threshold = NULL,
      parent = NULL,
      child = NULL
    )

#### Arguments

- `threshold`:

  minimum number of children, default `config$min_peptides_protein`

- `parent`:

  parent key column, default the first hierarchy key

- `child`:

  child key column, default the second hierarchy key

#### Returns

self (invisible)

------------------------------------------------------------------------

### `LFQData$remove_decoys()`

remove decoy / reverse-database entries (rows whose top-level hierarchy
id, e.g. protein_Id, is a decoy). Detection uses the configured
`pattern_decoys` unioned with the built-in defaults (see
[`is_decoy`](https://fgcz.github.io/prolfqua/reference/is_decoy.md));
decoys are always detectable even when no pattern is configured. Returns
a new decoy-free LFQData (self is not modified).

#### Usage

    LFQData$remove_decoys()

#### Returns

LFQData without decoys

------------------------------------------------------------------------

### `LFQData$remove_contaminants()`

remove contaminant entries (rows whose top-level hierarchy id, e.g.
protein_Id, is a contaminant). Detection uses the configured
`pattern_contaminants` unioned with the built-in defaults (see
[`is_contaminant`](https://fgcz.github.io/prolfqua/reference/is_contaminant.md)).
Returns a new contaminant-free LFQData (self is not modified).

#### Usage

    LFQData$remove_contaminants()

#### Returns

LFQData without contaminants

------------------------------------------------------------------------

### `LFQData$decoy_proportion()`

proportion of modelling-level keys (subject_id, e.g. protein or peptide
by `hierarchy_depth`) that are decoys - an empirical-FDR / target-decoy
QC signal. Decoy status derives from the top-level hierarchy id
(protein_Id).

#### Usage

    LFQData$decoy_proportion()

#### Returns

numeric in \[0, 1\]

------------------------------------------------------------------------

### `LFQData$complete_cases()`

some software is reporting NA's as 0, you must remove it from your data

#### Usage

    LFQData$complete_cases()

#### Returns

self

------------------------------------------------------------------------

### `LFQData$data_wide()`

converts the data to wide

#### Usage

    LFQData$data_wide(as.matrix = FALSE, value = NULL)

#### Arguments

- `as.matrix`:

  return as data.frame or matrix

- `value`:

  see possible lfqdata\$get_config()\$value_vars()

#### Returns

list with data, annotation, and configuration

------------------------------------------------------------------------

### `LFQData$factors()`

Annotation table

#### Usage

    LFQData$factors()

#### Returns

data.frame

------------------------------------------------------------------------

### `LFQData$hierarchy()`

Hierarchy table

#### Usage

    LFQData$hierarchy()

------------------------------------------------------------------------

### `LFQData$response()`

name of response variable

#### Usage

    LFQData$response()

#### Returns

data.frame

------------------------------------------------------------------------

### `LFQData$hierarchy_keys()`

return all hierarchy column names

#### Usage

    LFQData$hierarchy_keys()

------------------------------------------------------------------------

### `LFQData$relevant_hierarchy_keys()`

return hierarchy column names at current depth (alias for subject_id)

#### Usage

    LFQData$relevant_hierarchy_keys()

------------------------------------------------------------------------

### `LFQData$factor_keys()`

return all factor column names

#### Usage

    LFQData$factor_keys()

------------------------------------------------------------------------

### `LFQData$relevant_factor_keys()`

return factor column names at current depth

#### Usage

    LFQData$relevant_factor_keys()

------------------------------------------------------------------------

### `LFQData$sample_name()`

return sample name column

#### Usage

    LFQData$sample_name()

------------------------------------------------------------------------

### `LFQData$file_name()`

return file name column

#### Usage

    LFQData$file_name()

------------------------------------------------------------------------

### `LFQData$nr_children_col()`

return name of nr_children column

#### Usage

    LFQData$nr_children_col()

------------------------------------------------------------------------

### `LFQData$isotope_label()`

return isotope label column name

#### Usage

    LFQData$isotope_label()

------------------------------------------------------------------------

### `LFQData$data_long()`

return the tidy (long-format) data frame

#### Usage

    LFQData$data_long(na.omit = FALSE)

#### Arguments

- `na.omit`:

  if TRUE, remove rows with NA in response column

------------------------------------------------------------------------

### `LFQData$rename_response()`

new name of response variable. Renaming to the current response name
changes nothing; renaming to another existing column is an error.

#### Usage

    LFQData$rename_response(newname = "Intensity")

#### Arguments

- `newname`:

  default Intensity

#### Returns

self (invisible)

------------------------------------------------------------------------

### `LFQData$hierarchy_counts()`

number of elements at each level

#### Usage

    LFQData$hierarchy_counts()

------------------------------------------------------------------------

### `LFQData$summarize_hierarchy()`

e.g. number of peptides per protein etc

#### Usage

    LFQData$summarize_hierarchy()

#### Returns

data.frame

------------------------------------------------------------------------

### `LFQData$get_Plotter()`

get Plotter

#### Usage

    LFQData$get_Plotter()

#### Returns

LFQDataPlotter

------------------------------------------------------------------------

### `LFQData$get_Summariser()`

get Summariser

#### Usage

    LFQData$get_Summariser()

#### Returns

LFQDataSummarizer

------------------------------------------------------------------------

### `LFQData$get_Stats()`

Get
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md).
For more details see
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md).

#### Usage

    LFQData$get_Stats(stats = c("everything", "interaction", "all"))

#### Arguments

- `stats`:

  default interaction, computes statistics within interaction.

#### Returns

LFQDataStats

------------------------------------------------------------------------

### `LFQData$get_Transformer()`

get Stats

#### Usage

    LFQData$get_Transformer()

#### Returns

LFQDataTransformer

------------------------------------------------------------------------

### `LFQData$get_Aggregator()`

get Aggregator

#### Usage

    LFQData$get_Aggregator(method = "medpolish", ...)

#### Arguments

- `method`:

  aggregation method: "medpolish", "rlm", "topN", or "limpa"

- `...`:

  passed to aggregator constructor (e.g. prefix, N, func, dpc_slope)

#### Returns

AggregateMedpolish, AggregateRlm, AggregateTopN, or AggregateLimpa

------------------------------------------------------------------------

### `LFQData$clone()`

The objects of this class are cloneable with this method.

#### Usage

    LFQData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$filter_proteins_by_peptide_count()
#> Column added : nr_peptide_Id_IN_protein_Id
#> filter_proteins_by_peptide_count: 4 protein_Id removed (e.g. 7cbcrd~5725, 9VUkAq~4703, CGzoYe~2147, DoWup2~5896): fewer than 2 peptide_Id
#> See $provenance().
tmp <- lfqdata$data_wide()
stopifnot(nrow(tmp$data) == nrow(tmp$rowdata))
stopifnot(ncol(tmp$data) == nrow(tmp$annotation) + ncol(tmp$rowdata))

stopifnot("data.frame" %in% class(tmp$data))
tmp <- lfqdata$data_wide(as.matrix = TRUE)
stopifnot("matrix" %in% class(tmp$data))
stopifnot(lfqdata$is_transformed()==FALSE)
lfqdata$summarize_hierarchy()
#> # A tibble: 6 × 3
#>   protein_Id  isotopeLabel_n peptide_Id_n
#>   <chr>                <int>        <int>
#> 1 0EfVhX~0087              1            3
#> 2 BEJI92~5282              1            2
#> 3 Fl4JiV~8625              1            4
#> 4 HvIpHG~9079              1            2
#> 5 JcKVfU~9653              1            7
#> 6 SGIVBl~5782              1            6

lfqdata$response()
#> [1] "abundance"
lfqdata$rename_response("peptide.intensity")
lfqdata$response()
#> [1] "peptide.intensity"
stopifnot("LFQData" %in% class(lfqdata$get_copy()))
stopifnot("LFQDataTransformer" %in% class(lfqdata$get_Transformer()))
stopifnot("LFQDataStats" %in% class(lfqdata$get_Stats()))
stopifnot("LFQDataSummariser" %in% class(lfqdata$get_Summariser()))
stopifnot("LFQDataPlotter" %in% class(lfqdata$get_Plotter()))
stopifnot("AggregateMedpolish" %in% class(lfqdata$get_Aggregator("medpolish")))
#> Warning: You did not transform the intensities. medpolish works best with already variance stabilized intensities. Use LFQData$get_Transformer to transform the data: peptide.intensity

tmp <-lfqdata$get_sample(5, seed = 4)
stopifnot(nrow(tmp$hierarchy()) == 5)
```
