# Decorates LFQData with methods to aggregate protein intensities aggregates intensities

Decorates LFQData with methods to aggregate protein intensities
aggregates intensities

Decorates LFQData with methods to aggregate protein intensities
aggregates intensities

## See also

Other LFQData:
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

## Methods

### Public methods

- [`LFQDataAggregator$new()`](#method-LFQDataAggregator-new)

- [`LFQDataAggregator$medpolish()`](#method-LFQDataAggregator-medpolish)

- [`LFQDataAggregator$lmrob()`](#method-LFQDataAggregator-lmrob)

- [`LFQDataAggregator$mean_topN()`](#method-LFQDataAggregator-mean_topN)

- [`LFQDataAggregator$sum_topN()`](#method-LFQDataAggregator-sum_topN)

- [`LFQDataAggregator$plot()`](#method-LFQDataAggregator-plot)

- [`LFQDataAggregator$write_plots()`](#method-LFQDataAggregator-write_plots)

- [`LFQDataAggregator$clone()`](#method-LFQDataAggregator-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    LFQDataAggregator$new(lfq, prefix = "protein")

#### Arguments

- `lfq`:

  LFQData

- `prefix`:

  default protein

------------------------------------------------------------------------

### Method [`medpolish()`](https://rdrr.io/r/stats/medpolish.html)

aggregate using median polish

#### Usage

    LFQDataAggregator$medpolish()

#### Arguments

- `N`:

  top N by intensity

#### Returns

LFQData

------------------------------------------------------------------------

### Method `lmrob()`

aggregate using robust regression

#### Usage

    LFQDataAggregator$lmrob()

#### Arguments

- `N`:

  top N by intensity

#### Returns

LFQData

------------------------------------------------------------------------

### Method `mean_topN()`

aggregate topN using mean

#### Usage

    LFQDataAggregator$mean_topN(N = 3)

#### Arguments

- `N`:

  top N by intensity

#### Returns

LFQData

------------------------------------------------------------------------

### Method `sum_topN()`

aggregate topN using sum

#### Usage

    LFQDataAggregator$sum_topN(N = 3)

#### Arguments

- `N`:

  top N by intensity

#### Returns

LFQData

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

creates aggregation plots

#### Usage

    LFQDataAggregator$plot(subset = NULL, show.legend = FALSE)

#### Arguments

- `subset`:

  create plots for a subset of the data only, e.g. proteins with more
  then 2 peptides.

- `show.legend`:

  default FALSE

#### Returns

data.frame

------------------------------------------------------------------------

### Method `write_plots()`

writes plots to folder

#### Usage

    LFQDataAggregator$write_plots(
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

    LFQDataAggregator$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <-  prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
istar$config <- istar$config
data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
lfqdata <- LFQData$new(data, istar$config)

lfqTrans <- lfqdata$clone()$get_Transformer()
lfqTrans$log2()
#> Column added : log2_abundance
lfqTrans <- lfqTrans$robscale()$lfq
#> data is : TRUE
#> Warning: Expected 2 pieces. Additional pieces discarded in 336 rows [1, 2, 3, 4, 5, 6,
#> 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#> Joining with `by = join_by(sampleName, protein_Id, peptide_Id)`
lfqAggregator <- LFQDataAggregator$new(lfqTrans, "protein")

lfqAggregator$medpolish()
#> starting aggregation
pmed <- lfqAggregator$plot()
pmed$plots[[1]]
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).

lfqAggregator$lmrob()
#> starting aggregation
prob <- lfqAggregator$plot()
prob$plots[[1]]
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).


lfqCopy <- lfqdata$clone()
lfqCopy$is_transformed()
#> [1] FALSE
lfqAggregator <- LFQDataAggregator$new(lfqCopy, "protein")
lfqAggregator$sum_topN()
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank
pSum <- lfqAggregator$plot()
pSum$plots[[1]]
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).


lfqAggregator$mean_topN()
#> Joining with `by = join_by(protein_Id, peptide_Id)`
#> Columns added : srm_meanInt srm_meanIntRank
pMean <- lfqAggregator$plot()
pMean$plots[[1]]
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).

protPlotter <- lfqAggregator$lfq_agg$get_Plotter()
protPlotter$heatmap()

lfqAggregator$write_plots(tempdir())
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
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
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 8 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 6 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> Warning: Removed 3 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 9 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
```
