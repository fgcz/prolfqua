# Plot an LFQData object

Default plot: distribution of the response per sample (violin plot), on
a log scale when the response is not transformed. Other plots are
available from \`x\$get_Plotter()\`.

## Usage

``` r
# S3 method for class 'LFQData'
plot(x, ...)
```

## Arguments

- x:

  LFQData object

- ...:

  unused

## Value

ggplot

## See also

Other LFQData:
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md),
[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md),
[`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md),
[`LFQDataPlotter`](https://fgcz.github.io/prolfqua/reference/LFQDataPlotter.md),
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://fgcz.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md),
[`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md)

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)
lfqdata
#> LFQData: 12 samples
#> Response: abundance (not transformed) 
#> Hierarchy: protein_Id 10, peptide_Id 28 
#> Samples per group_: A 4, B 4, Ctrl 4
#> Response range: 13.8 to 48.5; missing 36 of 336 (10.7%)
#> Provenance: 0 failure and 0 other entries; see $provenance() 
plot(lfqdata)
#> Warning: Removed 36 rows containing non-finite outside the scale range
#> (`stat_ydensity()`).
#> Warning: Removed 36 rows containing non-finite outside the scale range
#> (`stat_summary()`).
```
