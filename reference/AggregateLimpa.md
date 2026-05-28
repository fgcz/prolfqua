# AggregateLimpa

AggregateLimpa

AggregateLimpa

## Value

An R6 class generator.

## Details

Aggregates peptide/precursor intensities to protein level using limpa's
probabilistic DPC-based quantification
([`limpa::dpcQuant`](https://rdrr.io/pkg/limpa/man/dpcQuant.html)). With
`impute_only = TRUE`, performs same-level imputation without aggregation
using
[`limpa::dpcQuantByRow`](https://rdrr.io/pkg/limpa/man/dpcQuant.html).

The output LFQData contains three value columns:

- intensity (`"limpa"`) — protein-level log-expression (complete, no
  NAs)

- standard error (in `config$opt_se`) — per-protein, per-sample
  posterior SEs

- observation count (in `config$nr_children`) — number of observed
  precursors

## See also

Other LFQData:
[`AggregateMedpolish`](https://wolski.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://wolski.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://wolski.github.io/prolfqua/reference/AggregateTopN.md),
[`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://wolski.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://wolski.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://wolski.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Public fields

- `lfq`:

  LFQData (deep cloned input)

- `lfq_agg`:

  aggregation result

- `prefix`:

  to use for aggregation results e.g. protein

- `dpc_result`:

  estimated DPC object from limpa::dpc

- `dpc_slope`:

  DPC slope parameter (default 0.8)

- `impute_only`:

  if TRUE use dpcQuantByRow (no aggregation)

## Methods

### Public methods

- [`AggregateLimpa$new()`](#method-AggregateLimpa-new)

- [`AggregateLimpa$aggregate()`](#method-AggregateLimpa-aggregate)

- [`AggregateLimpa$plot()`](#method-AggregateLimpa-plot)

- [`AggregateLimpa$write_plots()`](#method-AggregateLimpa-write_plots)

- [`AggregateLimpa$clone()`](#method-AggregateLimpa-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    AggregateLimpa$new(
      lfq,
      prefix = "protein",
      dpc_slope = 0.8,
      impute_only = FALSE
    )

#### Arguments

- `lfq`:

  LFQData with log2-transformed intensities

- `prefix`:

  default "protein"

- `dpc_slope`:

  DPC slope parameter passed to limpa (default 0.8)

- `impute_only`:

  if TRUE, use dpcQuantByRow instead of dpcQuant

------------------------------------------------------------------------

### Method [`aggregate()`](https://rdrr.io/r/stats/aggregate.html)

run limpa DPC-based aggregation (or imputation)

#### Usage

    AggregateLimpa$aggregate()

#### Returns

LFQData

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

creates aggregation plots (only for aggregation mode, not impute_only)

#### Usage

    AggregateLimpa$plot(subset = NULL, show.legend = FALSE)

#### Arguments

- `subset`:

  create plots for a subset of the data only

- `show.legend`:

  default FALSE

#### Returns

data.frame

------------------------------------------------------------------------

### Method `write_plots()`

writes plots to folder

#### Usage

    AggregateLimpa$write_plots(
      qcpath,
      subset = NULL,
      show.legend = FALSE,
      width = 6,
      height = 6
    )

#### Arguments

- `qcpath`:

  qcpath

- `subset`:

  write plots only for some

- `show.legend`:

  legend

- `width`:

  figure width

- `height`:

  figure height

#### Returns

file path

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    AggregateLimpa$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Column added : log2_abundance

agg <- AggregateLimpa$new(lfqdata, "protein")
agg$aggregate()
#> completing cases
agg$lfq_agg$data_wide()
#> $data
#> # A tibble: 10 × 14
#>    protein_Id  isotopeLabel  A_V1  A_V2  A_V3  A_V4  B_V1  B_V2  B_V3  B_V4
#>    <chr>       <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 0EfVhX~0087 light         4.40  4.32  4.26  4.34  4.74  4.68  4.73  4.85
#>  2 7cbcrd~5725 light         4.90  4.90  4.87  4.78  3.64  4.43  4.38  4.39
#>  3 9VUkAq~4703 light         4.10  4.12  4.25  4.24  4.04  3.59  4.13  4.14
#>  4 BEJI92~5282 light         4.31  4.37  4.36  4.37  4.54  4.47  4.47  4.40
#>  5 CGzoYe~2147 light         4.60  4.60  4.65  4.58  3.93  4.96  4.95  4.95
#>  6 DoWup2~5896 light         4.59  4.60  4.54  4.54  4.21  4.22  4.21  4.15
#>  7 Fl4JiV~8625 light         4.35  4.41  4.29  4.35  4.54  4.56  4.65  4.54
#>  8 HvIpHG~9079 light         4.15  4.29  4.10  4.24  4.74  4.78  4.78  4.74
#>  9 JcKVfU~9653 light         4.97  5.00  5.01  5.02  5.22  5.23  5.22  5.27
#> 10 SGIVBl~5782 light         4.67  4.61  4.76  4.73  4.73  4.69  4.71  4.74
#> # ℹ 4 more variables: Ctrl_V1 <dbl>, Ctrl_V2 <dbl>, Ctrl_V3 <dbl>,
#> #   Ctrl_V4 <dbl>
#> 
#> $annotation
#> # A tibble: 12 × 4
#>    sampleName sample  group_ isotopeLabel
#>    <chr>      <chr>   <chr>  <chr>       
#>  1 A_V1       A_V1    A      light       
#>  2 A_V2       A_V2    A      light       
#>  3 A_V3       A_V3    A      light       
#>  4 A_V4       A_V4    A      light       
#>  5 B_V1       B_V1    B      light       
#>  6 B_V2       B_V2    B      light       
#>  7 B_V3       B_V3    B      light       
#>  8 B_V4       B_V4    B      light       
#>  9 Ctrl_V1    Ctrl_V1 Ctrl   light       
#> 10 Ctrl_V2    Ctrl_V2 Ctrl   light       
#> 11 Ctrl_V3    Ctrl_V3 Ctrl   light       
#> 12 Ctrl_V4    Ctrl_V4 Ctrl   light       
#> 
#> $rowdata
#> # A tibble: 10 × 2
#>    protein_Id  isotopeLabel
#>    <chr>       <chr>       
#>  1 0EfVhX~0087 light       
#>  2 7cbcrd~5725 light       
#>  3 9VUkAq~4703 light       
#>  4 BEJI92~5282 light       
#>  5 CGzoYe~2147 light       
#>  6 DoWup2~5896 light       
#>  7 Fl4JiV~8625 light       
#>  8 HvIpHG~9079 light       
#>  9 JcKVfU~9653 light       
#> 10 SGIVBl~5782 light       
#> 
#> $config
#> <AnalysisConfiguration>
#>   Public:
#>     annotation_vars: function () 
#>     bin_resp: 
#>     clone: function (deep = FALSE) 
#>     factor_depth: 1
#>     factor_keys: function () 
#>     factor_keys_depth: function () 
#>     factors: list
#>     file_name: sample
#>     get_response: function () 
#>     hierarchy: list
#>     hierarchy_depth: 1
#>     hierarchy_keys: function (rev = FALSE) 
#>     hierarchy_keys_depth: function (names = TRUE) 
#>     id_required: function () 
#>     id_vars: function () 
#>     ident_q_value: qValue
#>     ident_score: 
#>     initialize: function () 
#>     is_response_transformed: TRUE
#>     isotope_label: isotopeLabel
#>     min_peptides_protein: 2
#>     norm_value: NULL
#>     nr_children: nr_children_protein_Id
#>     opt_mz: 
#>     opt_rt: 
#>     opt_se: limpa_se
#>     pop_response: function () 
#>     sample_name: sampleName
#>     sep: ~
#>     set_response: function (col_name) 
#>     value_vars: function () 
#>     work_intensity: limpa
#> 
```
