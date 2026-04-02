# compute pooled variance

following the documentation here:
https://online.stat.psu.edu/stat500/lesson/7/7.3/7.3.1/7.3.1.1

## Usage

``` r
pooled_V2(x)

pooled_V1(x)

compute_pooled(x, method = c("V1", "V2"))

poolvar(res1, config, method = c("V1", "V2"))
```

## Arguments

- x:

  data.frame

## Value

data.frame

## See also

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`lfq_power_t_test_quantiles_V2()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://wolski.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://wolski.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`summarize_stats()`](https://wolski.github.io/prolfqua/reference/summarize_stats.md)

Other stats:
[`INTERNAL_FUNCTIONS_BY_FAMILY`](https://wolski.github.io/prolfqua/reference/INTERNAL_FUNCTIONS_BY_FAMILY.md),
[`lfq_power_t_test_proteins()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_proteins.md),
[`lfq_power_t_test_quantiles()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles.md),
[`lfq_power_t_test_quantiles_V2()`](https://wolski.github.io/prolfqua/reference/lfq_power_t_test_quantiles_V2.md),
[`plot_stat_density()`](https://wolski.github.io/prolfqua/reference/plot_stat_density.md),
[`plot_stat_density_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_density_median.md),
[`plot_stat_violin()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin.md),
[`plot_stat_violin_median()`](https://wolski.github.io/prolfqua/reference/plot_stat_violin_median.md),
[`plot_stdv_vs_mean()`](https://wolski.github.io/prolfqua/reference/plot_stdv_vs_mean.md),
[`summarize_stats()`](https://wolski.github.io/prolfqua/reference/summarize_stats.md)

## Examples

``` r
x <- data.frame(nrMeasured =c(1,2,2), var = c(3,4,4), meanAbundance = c(3,3,3))
x <- data.frame(nrMeasured = c(1,2,1,1), var = c(NA, 0.0370, NA, NA),
  meanAbundance = c(-1.94,-1.46,-1.87,-1.45))
prolfqua:::pooled_V2(na.omit(x))
#>   n.groups n df        sd   var       sdT  mean
#> 1        1 2  1 0.1923538 0.037 0.1923538 -1.46
prolfqua:::pooled_V1(na.omit(x))
#>   n.groups n df        sd       sdT   var  mean
#> 1        1 2  1 0.1923538 0.1923538 0.037 -1.46
x <- x[1,, drop=FALSE]
x
#>   nrMeasured var meanAbundance
#> 1          1  NA         -1.94
na.omit(x)
#> [1] nrMeasured    var           meanAbundance
#> <0 rows> (or 0-length row.names)
prolfqua:::pooled_V2(na.omit(x))
#>   n.groups n df sd var sdT mean
#> 1        0 0  0  0   0 NaN  NaN
x <- data.frame(nrMeasured =c(1,2,2), var = c(3,4,4), meanAbundance = c(3,3,3))
x <- data.frame(nrMeasured = c(1,2,1,1), var = c(NA, 0.0370, NA, NA),
  meanAbundance = c(-1.94,-1.46,-1.87,-1.45))
compute_pooled(x)
#>   n.groups n df        sd       sdT   var  mean meanAll nrMeasured
#> 1        1 2  1 0.1923538 0.1923538 0.037 -1.46  -1.636          5
compute_pooled(x, method = "V2")
#>   n.groups n df        sd   var       sdT  mean meanAll nrMeasured
#> 1        1 2  1 0.1923538 0.037 0.1923538 -1.46  -1.636          5
y <- data.frame(dilution.=c("a","b","c"),
     nrReplicates = c(4,4,4), nrMeasured = c(0,0,1), sd =c(NA,NA,NA),
     var = c(NA,NA,NA),meanAbundance = c(NaN,NaN,NaN))
compute_pooled(y)
#>   n.groups n df  sd sdT var mean meanAll nrMeasured
#> 1        0 0  0 NaN NaN NaN  NaN     NaN          1
yb <- y |> dplyr::filter(nrMeasured > 1)

bb <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
config <- bb$config
data <- bb$data

res1 <- summarize_stats(data, config)
#> completing cases
pv <- poolvar(res1, config)
stopifnot(nrow(pv) == nrow(res1)/3)
```
