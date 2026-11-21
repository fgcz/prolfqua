# Aggregates e.g. protein abundances from peptide abundances

Aggregates e.g. protein abundances from peptide abundances

## Usage

``` r
estimate_intensity(lfqdata, method = c("medpolish", "rlm"), control = list())
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- method:

  "medpolish" (Tukey's median polish, column `medpolish`) or "rlm"
  (robust regression with
  [`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html), column `lmrob`)

- control:

  convergence controls merged over the defaults. medpolish: `maxiter`
  (default 10) and `eps` (default 0.01, the change of the sum of
  absolute residuals that stops the polish, on the log2 scale of the
  input); rlm: `maxit` (default 20) and `acc` (default 1e-4, relative
  change of the IRLS residuals).

## Value

returns list with data (data.frame), config (AnalysisConfiguration) and
failures (tibble with the subject `id` and `warning` text of every
warning raised while estimating one subject; the warnings are not
re-raised)

## See also

Other aggregation:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://fgcz.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`aggregate_intensity_top_n()`](https://fgcz.github.io/prolfqua/reference/aggregate_intensity_top_n.md),
[`plot_estimate()`](https://fgcz.github.io/prolfqua/reference/plot_estimate.md),
[`plot_hierarchies_line()`](https://fgcz.github.io/prolfqua/reference/plot_hierarchies_line.md)

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(dd$data, dd$config)
lfq <- lfq$get_Transformer()$log2()$lfq
#> Column added : log2_abundance
bbMed <- estimate_intensity(lfq, method = "medpolish")
bbRob <- estimate_intensity(lfq, method = "rlm")
nrow(bbMed$data)
#> [1] 116
nrow(bbRob$data)
#> [1] 116
xt <- dplyr::inner_join(bbMed$data, bbRob$data)
#> Joining with `by = join_by(protein_Id, sampleName, group_, sample,
#> isotopeLabel, nr_children_protein_Id)`
plot(xt$medpolish, xt$lmrob, log = "xy", pch = "*")
abline(0, 1, col = 2)

```
