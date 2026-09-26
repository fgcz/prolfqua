# AggregateRlm

AggregateRlm

AggregateRlm

## Value

An R6 class generator.

## Details

Aggregates peptide intensities to protein level using robust regression
(rlm). Works best with variance-stabilized (log-transformed)
intensities.

## See also

Other LFQData:
[`AggregateLimpa`](https://wolski.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://wolski.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateTopN`](https://wolski.github.io/prolfqua/reference/AggregateTopN.md),
[`AggregatorBase`](https://wolski.github.io/prolfqua/reference/AggregatorBase.md),
[`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://wolski.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://wolski.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://wolski.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Super class

[`prolfqua::AggregatorBase`](https://wolski.github.io/prolfqua/reference/AggregatorBase.md)
-\> `AggregateRlm`

## Methods

### Public methods

- [`AggregateRlm$aggregate()`](#method-AggregateRlm-aggregate)

- [`AggregateRlm$clone()`](#method-AggregateRlm-clone)

Inherited methods

- [`prolfqua::AggregatorBase$initialize()`](https://wolski.github.io/prolfqua/html/AggregatorBase.html#method-AggregatorBase-initialize)
- [`prolfqua::AggregatorBase$plot()`](https://wolski.github.io/prolfqua/html/AggregatorBase.html#method-AggregatorBase-plot)
- [`prolfqua::AggregatorBase$write_plots()`](https://wolski.github.io/prolfqua/html/AggregatorBase.html#method-AggregatorBase-write_plots)

------------------------------------------------------------------------

### Method [`aggregate()`](https://rdrr.io/r/stats/aggregate.html)

run robust regression aggregation

#### Usage

    AggregateRlm$aggregate()

#### Returns

LFQData

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
#> creating sampleName from file_name column
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
#> completing cases
p <- agg$plot()
p$plots[[1]]
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).

```
