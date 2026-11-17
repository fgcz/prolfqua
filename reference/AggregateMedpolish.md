# AggregateMedpolish

Aggregates peptide intensities to protein level using median polish.
Works best with variance-stabilized (log-transformed) intensities.

## Value

An R6 class generator.

## See also

Other LFQData:
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md),
[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md),
[`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://fgcz.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://fgcz.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md),
[`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md)

## Super class

[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md)
-\> `AggregateMedpolish`

## Methods

### Public methods

- [`AggregateMedpolish$aggregate()`](#method-AggregateMedpolish-aggregate)

- [`AggregateMedpolish$clone()`](#method-AggregateMedpolish-clone)

Inherited methods

- [`AggregatorBase$initialize()`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.html#method-initialize)
- [`AggregatorBase$plot()`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.html#method-plot)
- [`AggregatorBase$write_plots()`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.html#method-write_plots)

------------------------------------------------------------------------

### `AggregateMedpolish$aggregate()`

run median polish aggregation

#### Usage

    AggregateMedpolish$aggregate()

#### Returns

LFQData

------------------------------------------------------------------------

### `AggregateMedpolish$clone()`

The objects of this class are cloneable with this method.

#### Usage

    AggregateMedpolish$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
lfqdata <- LFQData$new(data, istar$config)
lfqTrans <- lfqdata$clone()$get_Transformer()$log2()$robscale()$lfq
#> Column added : log2_abundance

agg <- AggregateMedpolish$new(lfqTrans, "protein")
agg$aggregate()
p <- agg$plot()
p$plots[[1]]
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).

```
