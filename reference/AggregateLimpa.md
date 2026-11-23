# Aggregate or impute intensities with limpa's DPC quantification

Aggregates peptide/precursor intensities to protein level using limpa's
probabilistic DPC-based quantification
([`limpa::dpcQuant`](https://rdrr.io/pkg/limpa/man/dpcQuant.html)). With
`impute_only = TRUE`, performs same-level imputation without aggregation
using
[`limpa::dpcQuantByRow`](https://rdrr.io/pkg/limpa/man/dpcQuant.html).

The output LFQData contains three value columns:

- intensity (`"limpa"`) - protein-level log-expression (complete, no
  NAs)

- standard error (in `config$opt_se`) - per-protein, per-sample
  posterior SEs

- observation count (in `config$nr_children`) - number of observed
  precursors

limpa's iterative steps run with limpa's defaults:
[`limpa::dpc()`](https://rdrr.io/pkg/limpa/man/dpc.html) estimates the
detection probability curve in 2 iterations on a subset of 2000 rows,
and [`limpa::dpcQuant()`](https://rdrr.io/pkg/limpa/man/dpcQuant.html)
optimises each protein with its own settings. Only the slope bound
`dpc_slope` is exposed. limpa returns no convergence status for these
fits, so a changed setting could not be checked; the fit status of limpa
analyses reports `not_converged` as NA.

## Value

An R6 class generator.

## See also

Other LFQData:
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md),
[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md),
[`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://fgcz.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md),
[`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md),
[`lfqdata_to_summarized_experiment()`](https://fgcz.github.io/prolfqua/reference/lfqdata_to_summarized_experiment.md),
[`plot.LFQData()`](https://fgcz.github.io/prolfqua/reference/plot.LFQData.md)

## Super class

[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md)
-\> `AggregateLimpa`

## Public fields

- `dpc_result`:

  estimated DPC object from limpa::dpc

- `dpc_slope`:

  DPC slope parameter (default 0.8)

- `impute_only`:

  if TRUE use dpcQuantByRow (no aggregation)

## Methods

### Public methods

- [`AggregateLimpa$new()`](#method-AggregateLimpa-initialize)

- [`AggregateLimpa$aggregate()`](#method-AggregateLimpa-aggregate)

- [`AggregateLimpa$plot()`](#method-AggregateLimpa-plot)

- [`AggregateLimpa$write_plots()`](#method-AggregateLimpa-write_plots)

- [`AggregateLimpa$clone()`](#method-AggregateLimpa-clone)

------------------------------------------------------------------------

### `AggregateLimpa$new()`

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

### `AggregateLimpa$aggregate()`

run limpa DPC-based aggregation (or imputation)

#### Usage

    AggregateLimpa$aggregate()

#### Returns

LFQData

------------------------------------------------------------------------

### `AggregateLimpa$plot()`

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

### `AggregateLimpa$write_plots()`

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

### `AggregateLimpa$clone()`

The objects of this class are cloneable with this method.

#### Usage

    AggregateLimpa$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq
#> Transformed response column added: log2_abundance

agg <- AggregateLimpa$new(lfqdata, "protein")
agg$aggregate()
agg$lfq_agg$data_wide()
#> $data
#> # A tibble: 10 × 13
#>    protein_Id   A_V1  A_V2  A_V3  A_V4  B_V1  B_V2  B_V3  B_V4 Ctrl_V1 Ctrl_V2
#>    <chr>       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl>
#>  1 0EfVhX~0087  4.37  4.28  4.26  4.31  4.74  4.68  4.72  4.78    4.31    4.27
#>  2 7cbcrd~5725  4.90  4.90  4.86  4.77  3.55  4.43  4.38  4.39    4.24    3.95
#>  3 9VUkAq~4703  4.10  4.12  4.25  4.24  4.05  3.53  4.13  4.14    4.67    4.81
#>  4 BEJI92~5282  4.31  4.37  4.30  4.30  4.54  4.47  4.47  4.40    3.76    4.17
#>  5 CGzoYe~2147  4.60  4.60  4.65  4.58  3.66  4.96  4.95  4.95    4.82    4.98
#>  6 DoWup2~5896  4.59  4.60  4.54  4.54  4.21  4.22  4.21  4.16    4.27    4.27
#>  7 Fl4JiV~8625  4.32  4.37  4.29  4.27  4.54  4.56  4.60  4.51    4.40    4.36
#>  8 HvIpHG~9079  4.09  4.29  4.04  4.24  4.74  4.78  4.77  4.74    4.62    4.61
#>  9 JcKVfU~9653  4.97  5.00  5.01  5.02  5.22  5.23  5.22  5.27    5.14    5.09
#> 10 SGIVBl~5782  4.60  4.61  4.73  4.70  4.73  4.67  4.71  4.74    4.63    4.61
#> # ℹ 2 more variables: Ctrl_V3 <dbl>, Ctrl_V4 <dbl>
#> 
#> $annotation
#> # A tibble: 12 × 3
#>    sampleName sample  group_
#>    <chr>      <chr>   <chr> 
#>  1 A_V1       A_V1    A     
#>  2 A_V2       A_V2    A     
#>  3 A_V3       A_V3    A     
#>  4 A_V4       A_V4    A     
#>  5 B_V1       B_V1    B     
#>  6 B_V2       B_V2    B     
#>  7 B_V3       B_V3    B     
#>  8 B_V4       B_V4    B     
#>  9 Ctrl_V1    Ctrl_V1 Ctrl  
#> 10 Ctrl_V2    Ctrl_V2 Ctrl  
#> 11 Ctrl_V3    Ctrl_V3 Ctrl  
#> 12 Ctrl_V4    Ctrl_V4 Ctrl  
#> 
#> $rowdata
#> # A tibble: 10 × 1
#>    protein_Id 
#>    <chr>      
#>  1 0EfVhX~0087
#>  2 7cbcrd~5725
#>  3 9VUkAq~4703
#>  4 BEJI92~5282
#>  5 CGzoYe~2147
#>  6 DoWup2~5896
#>  7 Fl4JiV~8625
#>  8 HvIpHG~9079
#>  9 JcKVfU~9653
#> 10 SGIVBl~5782
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
#>     hierarchy_keys: function () 
#>     hierarchy_keys_depth: function (names = TRUE) 
#>     id_vars: function () 
#>     ident_q_value: qValue
#>     ident_score: 
#>     is_response_transformed: TRUE
#>     min_peptides_protein: 2
#>     norm_value: NULL
#>     nr_children: nr_children_protein_Id
#>     opt_mz: 
#>     opt_rt: 
#>     opt_se: limpa_se
#>     pattern_contaminants: NULL
#>     pattern_decoys: NULL
#>     pop_response: function () 
#>     sample_name: sampleName
#>     sep: ~
#>     set_response: function (col_name) 
#>     value_vars: function () 
#>     work_intensity: limpa
#> 
```
