# Aggregate peptide intensities by median polish

Aggregates peptide intensities to protein level using median polish.
Works best with variance-stabilized (log-transformed) intensities.

The aggregated data carry the peptide consistency of each protein in
each sample, from the median polish residuals of its peptides:
`consistency_sigma2` (the protein's residual variance, shrunk over
proteins with
[`limma::squeezeVar()`](https://rdrr.io/pkg/limma/man/squeezeVar.html)),
`consistency_factor` (the inverse of the sample's residual variance
relative to the protein's, shrunk over all protein and sample cells;
cells with one peptide get the prior) and `consistency_factor_rel` (the
factor divided by its mean over the protein's samples). They also carry
weight variants built from them (see
`AnalysisConfiguration$weight_columns`). The facades' default weights
multiply in `consistency_factor_rel` (see
`AnalysisConfiguration$use_peptide_consistency`).

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
[`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md),
[`lfqdata_to_summarized_experiment()`](https://fgcz.github.io/prolfqua/reference/lfqdata_to_summarized_experiment.md),
[`plot.LFQData()`](https://fgcz.github.io/prolfqua/reference/plot.LFQData.md)

## Super class

[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md)
-\> `AggregateMedpolish`

## Public fields

- `control`:

  convergence controls of
  [`stats::medpolish`](https://rdrr.io/r/stats/medpolish.html):
  `maxiter` (default 10) and `eps` (default 0.01; the polish stops when
  the sum of absolute residuals changes by less than `eps` times that
  sum)

- `consistency_priors`:

  prior df and variance of the peptide consistency, set by
  [`aggregate()`](https://rdrr.io/r/stats/aggregate.html): one row with
  `protein_df_prior`, `protein_var_prior` (protein level) and
  `cell_df_prior`, `cell_var_prior` (protein and sample level); NA when
  fewer than three proteins or cells have residual df

## Methods

### Public methods

- [`AggregateMedpolish$new()`](#method-AggregateMedpolish-initialize)

- [`AggregateMedpolish$aggregate()`](#method-AggregateMedpolish-aggregate)

- [`AggregateMedpolish$clone()`](#method-AggregateMedpolish-clone)

Inherited methods

- [`AggregatorBase$plot()`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.html#method-plot)
- [`AggregatorBase$write_plots()`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.html#method-write_plots)

------------------------------------------------------------------------

### `AggregateMedpolish$new()`

initialize

#### Usage

    AggregateMedpolish$new(lfq, prefix = "protein", control = list())

#### Arguments

- `lfq`:

  LFQData

- `prefix`:

  default protein

- `control`:

  named list merged over the defaults; convergence controls of
  [`stats::medpolish`](https://rdrr.io/r/stats/medpolish.html):
  `maxiter` (default 10) and `eps` (default 0.01; the polish stops when
  the sum of absolute residuals changes by less than `eps` times that
  sum)

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
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 100))
lfqdata <- LFQData$new(data, istar$config)
lfqTrans <- lfqdata$clone()$get_Transformer()$log2()$robscale()$lfq
#> Transformed response column added: log2_abundance

agg <- AggregateMedpolish$new(lfqTrans, "protein")
agg$aggregate()
p <- agg$plot()
p$plots[[1]]

```
