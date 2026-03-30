# AggregateRlm

AggregateRlm

AggregateRlm

## Details

Aggregates peptide intensities to protein level using robust regression
(rlm). Works best with variance-stabilized (log-transformed)
intensities.

## See also

Other LFQData:
[`AggregateLimpa`](https://wolski.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://wolski.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateTopN`](https://wolski.github.io/prolfqua/reference/AggregateTopN.md),
[`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://wolski.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://wolski.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://wolski.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Public fields

- `lfq`:

  LFQData

- `lfq_agg`:

  aggregation result

- `prefix`:

  to use for aggregation results e.g. protein

## Methods

### Public methods

- [`AggregateRlm$new()`](#method-AggregateRlm-new)

- [`AggregateRlm$aggregate()`](#method-AggregateRlm-aggregate)

- [`AggregateRlm$plot()`](#method-AggregateRlm-plot)

- [`AggregateRlm$write_plots()`](#method-AggregateRlm-write_plots)

- [`AggregateRlm$clone()`](#method-AggregateRlm-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    AggregateRlm$new(lfq, prefix = "protein")

#### Arguments

- `lfq`:

  LFQData

- `prefix`:

  default protein

------------------------------------------------------------------------

### Method [`aggregate()`](https://rdrr.io/r/stats/aggregate.html)

run robust regression aggregation

#### Usage

    AggregateRlm$aggregate()

#### Returns

LFQData

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

creates aggregation plots

#### Usage

    AggregateRlm$plot(subset = NULL, show.legend = FALSE)

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

    AggregateRlm$write_plots(
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

    AggregateRlm$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
lfqdata <- LFQData$new(data, istar$config)
lfqTrans <- lfqdata$clone()$get_Transformer()$log2()$robscale()$lfq
#> Column added : log2_abundance
#> data is : TRUE
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id, peptide_Id)`

agg <- AggregateRlm$new(lfqTrans, "protein")
agg$aggregate()
#> starting aggregation
p <- agg$plot()
p$plots[[1]]
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).

```
