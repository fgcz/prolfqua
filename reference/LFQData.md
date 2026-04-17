# LFQData R6 class

LFQData R6 class

LFQData R6 class

## Missing Data Assumptions

The filtering and imputation methods in this package assume that missing
values are Missing Completely At Random (MCAR) or Missing At Random
(MAR). Abundance-dependent missingness (MNAR), which is common in DDA
proteomics, is not modelled. Users should be aware that MNAR can bias
fold-change estimates and inflate false discovery rates.

## See also

Other LFQData:
[`AggregateLimpa`](https://wolski.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://wolski.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://wolski.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://wolski.github.io/prolfqua/reference/AggregateTopN.md),
[`LFQDataPlotter`](https://wolski.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://wolski.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://wolski.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Public fields

- `prefix`:

  e.g. "peptide\_", "protein\_", "compound\_"

## Methods

### Public methods

- [`LFQData$new()`](#method-LFQData-new)

- [`LFQData$set_data()`](#method-LFQData-set_data)

- [`LFQData$get_config()`](#method-LFQData-get_config)

- [`LFQData$set_config_value()`](#method-LFQData-set_config_value)

- [`LFQData$get_copy()`](#method-LFQData-get_copy)

- [`LFQData$get_sample()`](#method-LFQData-get_sample)

- [`LFQData$get_subset()`](#method-LFQData-get_subset)

- [`LFQData$subject_id()`](#method-LFQData-subject_id)

- [`LFQData$is_transformed()`](#method-LFQData-is_transformed)

- [`LFQData$remove_small_intensities()`](#method-LFQData-remove_small_intensities)

- [`LFQData$filter_proteins_by_peptide_count()`](#method-LFQData-filter_proteins_by_peptide_count)

- [`LFQData$omit_na()`](#method-LFQData-omit_na)

- [`LFQData$complete_cases()`](#method-LFQData-complete_cases)

- [`LFQData$data_wide()`](#method-LFQData-data_wide)

- [`LFQData$to_wide()`](#method-LFQData-to_wide)

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

- [`LFQData$get_data()`](#method-LFQData-get_data)

- [`LFQData$rename_response()`](#method-LFQData-rename_response)

- [`LFQData$hierarchy_counts()`](#method-LFQData-hierarchy_counts)

- [`LFQData$summarize_hierarchy()`](#method-LFQData-summarize_hierarchy)

- [`LFQData$get_Plotter()`](#method-LFQData-get_Plotter)

- [`LFQData$get_Summariser()`](#method-LFQData-get_Summariser)

- [`LFQData$get_Stats()`](#method-LFQData-get_Stats)

- [`LFQData$get_Transformer()`](#method-LFQData-get_Transformer)

- [`LFQData$get_Aggregator()`](#method-LFQData-get_Aggregator)

- [`LFQData$filter_difference()`](#method-LFQData-filter_difference)

- [`LFQData$clone()`](#method-LFQData-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    LFQData$new(data, config, prefix = "ms_", setup = FALSE)

#### Arguments

- `data`:

  data.frame

- `config`:

  configuration

- `prefix`:

  will be use as output prefix

- `setup`:

  is data setup needed, default = FALSE, if TRUE, calls
  [`setup_analysis`](https://wolski.github.io/prolfqua/reference/setup_analysis.md)
  on data first.

- `is_pep`:

  todo

------------------------------------------------------------------------

### Method `set_data()`

set data (replaces the internal data frame)

#### Usage

    LFQData$set_data(new_data)

#### Arguments

- `new_data`:

  data.frame

#### Returns

self (invisible)

------------------------------------------------------------------------

### Method `get_config()`

return the AnalysisConfiguration object

#### Usage

    LFQData$get_config()

#### Returns

AnalysisConfiguration

------------------------------------------------------------------------

### Method `set_config_value()`

set a config field value

#### Usage

    LFQData$set_config_value(field, value)

#### Arguments

- `field`:

  character — field name

- `value`:

  the value to set

------------------------------------------------------------------------

### Method `get_copy()`

get deep copy

#### Usage

    LFQData$get_copy()

------------------------------------------------------------------------

### Method `get_sample()`

samples subset of data

#### Usage

    LFQData$get_sample(size = 100, seed = NULL)

#### Arguments

- `size`:

  size of subset default 100

- `seed`:

  set seed

------------------------------------------------------------------------

### Method `get_subset()`

get subset of data

#### Usage

    LFQData$get_subset(x)

#### Arguments

- `x`:

  data frame with columns containing subject_id

------------------------------------------------------------------------

### Method `subject_id()`

get subject ID columns

#### Usage

    LFQData$subject_id()

------------------------------------------------------------------------

### Method `is_transformed()`

is data transformed

#### Usage

    LFQData$is_transformed(is_transformed)

#### Arguments

- `is_transformed`:

  logical

#### Returns

logical

------------------------------------------------------------------------

### Method [`remove_small_intensities()`](https://wolski.github.io/prolfqua/reference/remove_small_intensities.md)

some software is reporting NA's as 0, you must remove it from your data

#### Usage

    LFQData$remove_small_intensities(threshold = 4)

#### Arguments

- `threshold`:

  default 4.

#### Returns

self

------------------------------------------------------------------------

### Method [`filter_proteins_by_peptide_count()`](https://wolski.github.io/prolfqua/reference/filter_proteins_by_peptide_count.md)

remove proteins with less than X peptides

#### Usage

    LFQData$filter_proteins_by_peptide_count()

#### Returns

self

------------------------------------------------------------------------

### Method `omit_na()`

Omit NA from intensities per hierarchy (e.g. protein or peptide), idea
is to use it for normalization For instance if a peptide has a missing
value in more then nrNA of the samples within a condition it will be
removed

#### Usage

    LFQData$omit_na(nr_na = 0, factor_depth = NULL)

#### Arguments

- `nr_na`:

  number of NA values

- `factor_depth`:

  control whether \`nr_na\` is applied per condition or more globally,
  e.g. \`factor_depth = 0\` means per experiment

#### Returns

LFQData with NA omitted.

------------------------------------------------------------------------

### Method [`complete_cases()`](https://wolski.github.io/prolfqua/reference/complete_cases.md)

some software is reporting NA's as 0, you must remove it from your data

#### Usage

    LFQData$complete_cases()

#### Arguments

- `threshold`:

  default 4.

#### Returns

self

------------------------------------------------------------------------

### Method `data_wide()`

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

### Method `to_wide()`

deprecated — use data_wide() instead

#### Usage

    LFQData$to_wide(as.matrix = FALSE, value = NULL)

#### Arguments

- `as.matrix`:

  return as data.frame or matrix

- `value`:

  see possible lfqdata\$get_config()\$value_vars()

------------------------------------------------------------------------

### Method `factors()`

Annotation table

#### Usage

    LFQData$factors()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `hierarchy()`

Hierarchy table

#### Usage

    LFQData$hierarchy()

------------------------------------------------------------------------

### Method `response()`

name of response variable

#### Usage

    LFQData$response()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `hierarchy_keys()`

return all hierarchy column names

#### Usage

    LFQData$hierarchy_keys()

------------------------------------------------------------------------

### Method `relevant_hierarchy_keys()`

return hierarchy column names at current depth (alias for subject_id)

#### Usage

    LFQData$relevant_hierarchy_keys()

------------------------------------------------------------------------

### Method `factor_keys()`

return all factor column names

#### Usage

    LFQData$factor_keys()

------------------------------------------------------------------------

### Method `relevant_factor_keys()`

return factor column names at current depth

#### Usage

    LFQData$relevant_factor_keys()

------------------------------------------------------------------------

### Method `sample_name()`

return sample name column

#### Usage

    LFQData$sample_name()

------------------------------------------------------------------------

### Method `file_name()`

return file name column

#### Usage

    LFQData$file_name()

------------------------------------------------------------------------

### Method `nr_children_col()`

return name of nr_children column

#### Usage

    LFQData$nr_children_col()

------------------------------------------------------------------------

### Method `isotope_label()`

return isotope label column name

#### Usage

    LFQData$isotope_label()

------------------------------------------------------------------------

### Method `data_long()`

return the tidy (long-format) data frame

#### Usage

    LFQData$data_long(na.omit = FALSE)

#### Arguments

- `na.omit`:

  if TRUE, remove rows with NA in response column

------------------------------------------------------------------------

### Method `get_data()`

deprecated — use data_long() instead

#### Usage

    LFQData$get_data()

------------------------------------------------------------------------

### Method `rename_response()`

new name of response variable

#### Usage

    LFQData$rename_response(newname = "Intensity")

#### Arguments

- `newname`:

  default Intensity

------------------------------------------------------------------------

### Method [`hierarchy_counts()`](https://wolski.github.io/prolfqua/reference/hierarchy_counts.md)

number of elements at each level

#### Usage

    LFQData$hierarchy_counts()

------------------------------------------------------------------------

### Method [`summarize_hierarchy()`](https://wolski.github.io/prolfqua/reference/summarize_hierarchy.md)

e.g. number of peptides per protein etc

#### Usage

    LFQData$summarize_hierarchy()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `get_Plotter()`

get Plotter

#### Usage

    LFQData$get_Plotter()

#### Returns

LFQDataPlotter

------------------------------------------------------------------------

### Method `get_Summariser()`

get Summariser

#### Usage

    LFQData$get_Summariser()

#### Returns

LFQDataSummarizer

------------------------------------------------------------------------

### Method `get_Stats()`

Get
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md).
For more details see
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md).

#### Usage

    LFQData$get_Stats(stats = c("everything", "interaction", "all"))

#### Arguments

- `stats`:

  default interaction, computes statistics within interaction.

#### Returns

LFQDataStats

------------------------------------------------------------------------

### Method `get_Transformer()`

get Stats

#### Usage

    LFQData$get_Transformer()

#### Returns

LFQDataTransformer

------------------------------------------------------------------------

### Method `get_Aggregator()`

get Aggregator

#### Usage

    LFQData$get_Aggregator(method = "medpolish", ...)

#### Arguments

- `method`:

  aggregation method: "medpolish", "rlm", or "topN"

- `...`:

  passed to aggregator constructor (e.g. prefix, N, func)

#### Returns

AggregateMedpolish, AggregateRlm, or AggregateTopN

------------------------------------------------------------------------

### Method [`filter_difference()`](https://wolski.github.io/prolfqua/reference/filter_difference.md)

get difference of self with other if other is subset of self

#### Usage

    LFQData$filter_difference(other)

#### Arguments

- `other`:

  a filtered LFQData set

#### Details

Use to compare filtering results obtained from self, e.g. which proteins
and peptides were removed (other)

#### Returns

LFQData

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    LFQData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata$filter_proteins_by_peptide_count()
#> removing proteins with less than: 2 peptpides
#> Column added : nr_peptide_Id_IN_protein_Id
tmp <- lfqdata$data_wide()
testthat::expect_equal(nrow(tmp$data) , nrow(tmp$rowdata))
testthat::expect_equal(ncol(tmp$data) , nrow(tmp$annotation) + ncol(tmp$rowdata))

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

# filter for missing values

f1 <- lfqdata$omit_na(nr_na = 0)
#> Joining with `by = join_by(protein_Id, peptide_Id)`
stopifnot(f1$hierarchy_counts() <= lfqdata$hierarchy_counts())

f2 <- lfqdata$omit_na(factor_depth = 0)
#> Joining with `by = join_by(protein_Id, peptide_Id)`
stopifnot(f2$hierarchy_counts() <= lfqdata$hierarchy_counts())

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

lfqdata2 <- lfqdata$get_copy()
lfqdata2$set_data(lfqdata2$data_long()[1:100, ])
res <- lfqdata$filter_difference(lfqdata2)
stopifnot(nrow(res$data_long()) == nrow(lfqdata$data_long()) - 100)

tmp <- lfqdata$get_sample(5, seed = 4)
#> Sampling 5protein_Id
#> Joining with `by = join_by(protein_Id)`
stopifnot(nrow(tmp$hierarchy()) == 5)
```
