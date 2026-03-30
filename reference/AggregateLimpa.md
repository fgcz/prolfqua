# AggregateLimpa

AggregateLimpa

AggregateLimpa

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
if (FALSE) { # \dontrun{
istar <- prolfqua::sim_lfq_data_peptide_config()
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata <- lfqdata$get_Transformer()$log2()$lfq

agg <- AggregateLimpa$new(lfqdata, "protein")
agg$aggregate()
agg$lfq_agg$to_wide()
} # }
```
