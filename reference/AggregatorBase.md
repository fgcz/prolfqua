# Base class of the peptide to protein aggregators

Holds the fields, construction, and plotting shared by
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md)
and
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md).
Subclasses implement
[`aggregate()`](https://rdrr.io/r/stats/aggregate.html).

## Value

An R6 class generator.

## See also

Other LFQData:
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md),
[`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://fgcz.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://fgcz.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Public fields

- `lfq`:

  LFQData

- `lfq_agg`:

  aggregation result

- `prefix`:

  to use for aggregation results e.g. protein

## Methods

### Public methods

- [`AggregatorBase$new()`](#method-AggregatorBase-initialize)

- [`AggregatorBase$plot()`](#method-AggregatorBase-plot)

- [`AggregatorBase$write_plots()`](#method-AggregatorBase-write_plots)

- [`AggregatorBase$clone()`](#method-AggregatorBase-clone)

------------------------------------------------------------------------

### `AggregatorBase$new()`

initialize

#### Usage

    AggregatorBase$new(lfq, prefix = "protein")

#### Arguments

- `lfq`:

  LFQData

- `prefix`:

  default protein

------------------------------------------------------------------------

### `AggregatorBase$plot()`

creates aggregation plots

#### Usage

    AggregatorBase$plot(subset = NULL, show.legend = FALSE)

#### Arguments

- `subset`:

  create plots for a subset of the data only

- `show.legend`:

  default FALSE

#### Returns

data.frame

------------------------------------------------------------------------

### `AggregatorBase$write_plots()`

writes plots to folder

#### Usage

    AggregatorBase$write_plots(
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

### `AggregatorBase$clone()`

The objects of this class are cloneable with this method.

#### Usage

    AggregatorBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
