# AggregateTopN

AggregateTopN

AggregateTopN

## Details

Aggregates peptide intensities to protein level using top N peptides.
Works with raw (untransformed) intensities.

## See also

Other LFQData:
[`AggregateMedpolish`](https://wolski.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://wolski.github.io/prolfqua/reference/AggregateRlm.md),
[`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md),
[`LFQDataImp`](https://wolski.github.io/prolfqua/reference/LFQDataImp.md),
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

- `N`:

  top N peptides by intensity

- `func`:

  aggregation function name: "sum" or "mean"

## Methods

### Public methods

- [`AggregateTopN$new()`](#method-AggregateTopN-new)

- [`AggregateTopN$aggregate()`](#method-AggregateTopN-aggregate)

- [`AggregateTopN$plot()`](#method-AggregateTopN-plot)

- [`AggregateTopN$write_plots()`](#method-AggregateTopN-write_plots)

- [`AggregateTopN$clone()`](#method-AggregateTopN-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    AggregateTopN$new(lfq, prefix = "protein", N = 3, func = "sum")

#### Arguments

- `lfq`:

  LFQData

- `prefix`:

  default protein

- `N`:

  top N peptides (default 3)

- `func`:

  "sum" or "mean" (default "sum")

------------------------------------------------------------------------

### Method [`aggregate()`](https://rdrr.io/r/stats/aggregate.html)

run top N aggregation

#### Usage

    AggregateTopN$aggregate()

#### Returns

LFQData

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

creates aggregation plots

#### Usage

    AggregateTopN$plot(subset = NULL, show.legend = FALSE)

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

    AggregateTopN$write_plots(
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

    AggregateTopN$clone(deep = FALSE)

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

agg <- AggregateTopN$new(lfqdata, "protein", N = 3, func = "sum")
agg$aggregate()
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank
p <- agg$plot()
p$plots[[1]]
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).


agg_mean <- AggregateTopN$new(lfqdata, "protein", N = 3, func = "mean")
agg_mean$aggregate()
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank
protPlotter <- agg_mean$lfq_agg$get_Plotter()
protPlotter$heatmap()

agg_mean$write_plots(tempdir())
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 5 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 8 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_line()`).
```
