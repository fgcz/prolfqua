# Aggregate peptide intensities by robust linear regression

Aggregates peptide intensities to protein level using robust regression
(rlm). Works best with variance-stabilized (log-transformed)
intensities.

The aggregated data carry the peptide consistency columns and weight
variants of
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
computed from the rlm residuals.

## Value

An R6 class generator.

## See also

Other LFQData:
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
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
-\> `AggregateRlm`

## Public fields

- `control`:

  convergence controls of
  [`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html): `maxit`
  (default 20) and `acc` (default 1e-4, relative change of the IRLS
  residuals)

- `consistency_priors`:

  prior df and variance of the peptide consistency, set by
  [`aggregate()`](https://rdrr.io/r/stats/aggregate.html); see
  [`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md)

## Methods

### Public methods

- [`AggregateRlm$new()`](#method-AggregateRlm-initialize)

- [`AggregateRlm$aggregate()`](#method-AggregateRlm-aggregate)

- [`AggregateRlm$clone()`](#method-AggregateRlm-clone)

Inherited methods

- [`AggregatorBase$plot()`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.html#method-plot)
- [`AggregatorBase$write_plots()`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.html#method-write_plots)

------------------------------------------------------------------------

### `AggregateRlm$new()`

initialize

#### Usage

    AggregateRlm$new(lfq, prefix = "protein", control = list())

#### Arguments

- `lfq`:

  LFQData

- `prefix`:

  default protein

- `control`:

  named list merged over the defaults; convergence controls of
  [`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html): `maxit`
  (default 20) and `acc` (default 1e-4, relative change of the IRLS
  residuals)

------------------------------------------------------------------------

### `AggregateRlm$aggregate()`

run robust regression aggregation

#### Usage

    AggregateRlm$aggregate()

#### Returns

LFQData

------------------------------------------------------------------------

### `AggregateRlm$clone()`

The objects of this class are cloneable with this method.

#### Usage

    AggregateRlm$clone(deep = FALSE)

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
data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
lfqdata <- LFQData$new(data, istar$config)
lfqTrans <- lfqdata$clone()$get_Transformer()$log2()$robscale()$lfq
#> Transformed response column added: log2_abundance

agg <- AggregateRlm$new(lfqTrans, "protein")
agg$aggregate()
p <- agg$plot()
p$plots[[1]]

```
