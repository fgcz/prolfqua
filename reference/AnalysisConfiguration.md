# Analysis Configuration

Analysis Configuration — holds all table annotations, hierarchy
definitions, factor definitions, and analysis parameters in a single
flat object.

## See also

Other configuration:
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

- `sep`:

  separator to use when uniting columns is necessary

- `fileName`:

  column name of column containing raw file names

- `sampleName`:

  (will be generated from factors or fileName)

- `normValue`:

  optional column with normalization values (e.g., Creatinine)

- `isotopeLabel`:

  which column contains the isotope label (e.g. heavy or light), or
  light only if LFQ.

- `ident_qValue`:

  column name with identification QValues (smaller better)

- `ident_Score`:

  column with identification score (larger better)

- `opt_rt`:

  optional column with rt information

- `opt_mz`:

  optional column with mz information

- `nr_children`:

  optional column containing for instance the number of peptides

- `workIntensity`:

  column which contains the intensities

- `is_response_transformed`:

  are the intensities transformed for constant variance

- `bin_resp`:

  column with encoded missing information

- `factors`:

  Names of columns containing factors (annotations)

- `factorDepth`:

  number of relevant factors (used by plotting functions etc)

- `hierarchy`:

  list with columns describing the measurement hierarchy (i.e. protein
  peptide precursor fragment)

- `hierarchyDepth`:

  At which depth do you want to model i.e. protein than 1

- `min_peptides_protein`:

  minimum number of peptides per protein

## Active bindings

- `table`:

  deprecated, use config directly. Returns self for backwards
  compatibility.

- `parameter`:

  deprecated, use config directly. Returns self for backwards
  compatibility.

## Methods

### Public methods

- [`AnalysisConfiguration$new()`](#method-AnalysisConfiguration-new)

- [`AnalysisConfiguration$set_response()`](#method-AnalysisConfiguration-set_response)

- [`AnalysisConfiguration$get_response()`](#method-AnalysisConfiguration-get_response)

- [`AnalysisConfiguration$pop_response()`](#method-AnalysisConfiguration-pop_response)

- [`AnalysisConfiguration$factor_keys()`](#method-AnalysisConfiguration-factor_keys)

- [`AnalysisConfiguration$factor_keys_depth()`](#method-AnalysisConfiguration-factor_keys_depth)

- [`AnalysisConfiguration$hierarchy_keys()`](#method-AnalysisConfiguration-hierarchy_keys)

- [`AnalysisConfiguration$hierarchyKeys()`](#method-AnalysisConfiguration-hierarchyKeys)

- [`AnalysisConfiguration$hierarchy_keys_depth()`](#method-AnalysisConfiguration-hierarchy_keys_depth)

- [`AnalysisConfiguration$hkeysDepth()`](#method-AnalysisConfiguration-hkeysDepth)

- [`AnalysisConfiguration$id_required()`](#method-AnalysisConfiguration-id_required)

- [`AnalysisConfiguration$id_vars()`](#method-AnalysisConfiguration-id_vars)

- [`AnalysisConfiguration$value_vars()`](#method-AnalysisConfiguration-value_vars)

- [`AnalysisConfiguration$annotation_vars()`](#method-AnalysisConfiguration-annotation_vars)

- [`AnalysisConfiguration$clone()`](#method-AnalysisConfiguration-clone)

------------------------------------------------------------------------

### Method `new()`

create AnalysisConfiguration

#### Usage

    AnalysisConfiguration$new(
      analysisTableAnnotation = NULL,
      analysisParameter = NULL
    )

#### Arguments

- `analysisTableAnnotation`:

  optional, for backwards compatibility with old constructor

- `analysisParameter`:

  optional, for backwards compatibility with old constructor

------------------------------------------------------------------------

### Method `set_response()`

Add name of intensity column

#### Usage

    AnalysisConfiguration$set_response(colName)

#### Arguments

- `colName`:

  name of intensity column

------------------------------------------------------------------------

### Method `get_response()`

Get name of working intensity column

#### Usage

    AnalysisConfiguration$get_response()

------------------------------------------------------------------------

### Method `pop_response()`

Remove last name in array of working intensity column names

#### Usage

    AnalysisConfiguration$pop_response()

------------------------------------------------------------------------

### Method `factor_keys()`

Get factor keys

#### Usage

    AnalysisConfiguration$factor_keys()

#### Returns

array with keys

------------------------------------------------------------------------

### Method `factor_keys_depth()`

Get factor keys till factorDepth

#### Usage

    AnalysisConfiguration$factor_keys_depth()

------------------------------------------------------------------------

### Method `hierarchy_keys()`

get hierarchy keys

#### Usage

    AnalysisConfiguration$hierarchy_keys(rev = FALSE)

#### Arguments

- `rev`:

  return in reverse order

#### Returns

array of column names

------------------------------------------------------------------------

### Method `hierarchyKeys()`

get hierarchy keys (deprecated alias for hierarchy_keys)

#### Usage

    AnalysisConfiguration$hierarchyKeys(rev = FALSE)

#### Arguments

- `rev`:

  return in reverse order

#### Returns

array of column names

------------------------------------------------------------------------

### Method `hierarchy_keys_depth()`

get hierarchy keys up to depth

#### Usage

    AnalysisConfiguration$hierarchy_keys_depth(names = TRUE)

#### Arguments

- `names`:

  if TRUE names only if FALSE key value pairs

#### Returns

array of column names

------------------------------------------------------------------------

### Method `hkeysDepth()`

get hierarchy keys up to depth (deprecated alias for
hierarchy_keys_depth)

#### Usage

    AnalysisConfiguration$hkeysDepth(names = TRUE)

#### Arguments

- `names`:

  if TRUE names only if FALSE key value pairs

#### Returns

array of column names

------------------------------------------------------------------------

### Method `id_required()`

Id Columns which must be in the input data frame

#### Usage

    AnalysisConfiguration$id_required()

#### Returns

character array

------------------------------------------------------------------------

### Method `id_vars()`

get names of columns annotating values (e.g. intensities)

#### Usage

    AnalysisConfiguration$id_vars()

#### Returns

character array

------------------------------------------------------------------------

### Method `value_vars()`

get names of columns containing observations e.g. (intensity, qValue, mz
or rt)

#### Usage

    AnalysisConfiguration$value_vars()

------------------------------------------------------------------------

### Method `annotation_vars()`

get names of columns with sample annotations

#### Usage

    AnalysisConfiguration$annotation_vars()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    AnalysisConfiguration$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- istar$config
stopifnot("AnalysisConfiguration" %in% class(config))
stopifnot(length(config$hierarchy_keys()) > 0)
stopifnot(length(config$factor_keys()) > 0)
stopifnot(length(config$get_response()) == 1)
```
