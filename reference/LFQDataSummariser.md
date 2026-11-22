# Summarize LFQData

Summarize LFQData

## Value

An R6 class generator.

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
[`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md),
[`lfqdata_to_summarized_experiment()`](https://fgcz.github.io/prolfqua/reference/lfqdata_to_summarized_experiment.md),
[`plot.LFQData()`](https://fgcz.github.io/prolfqua/reference/plot.LFQData.md)

## Public fields

- `lfq`:

  LFQData

## Methods

### Public methods

- [`LFQDataSummariser$new()`](#method-LFQDataSummariser-initialize)

- [`LFQDataSummariser$hierarchy_counts()`](#method-LFQDataSummariser-hierarchy_counts)

- [`LFQDataSummariser$hierarchy_counts_sample()`](#method-LFQDataSummariser-hierarchy_counts_sample)

- [`LFQDataSummariser$plot_hierarchy_counts_sample()`](#method-LFQDataSummariser-plot_hierarchy_counts_sample)

- [`LFQDataSummariser$interaction_missing_stats()`](#method-LFQDataSummariser-interaction_missing_stats)

- [`LFQDataSummariser$upset_interaction_missing_stats()`](#method-LFQDataSummariser-upset_interaction_missing_stats)

- [`LFQDataSummariser$missingness_per_group()`](#method-LFQDataSummariser-missingness_per_group)

- [`LFQDataSummariser$missingness_per_group_cumsum()`](#method-LFQDataSummariser-missingness_per_group_cumsum)

- [`LFQDataSummariser$plot_missingness_per_group()`](#method-LFQDataSummariser-plot_missingness_per_group)

- [`LFQDataSummariser$plot_missingness_per_group_cumsum()`](#method-LFQDataSummariser-plot_missingness_per_group_cumsum)

- [`LFQDataSummariser$percentage_abundance()`](#method-LFQDataSummariser-percentage_abundance)

- [`LFQDataSummariser$clone()`](#method-LFQDataSummariser-clone)

------------------------------------------------------------------------

### `LFQDataSummariser$new()`

initialize

#### Usage

    LFQDataSummariser$new(lfqdata)

#### Arguments

- `lfqdata`:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

------------------------------------------------------------------------

### `LFQDataSummariser$hierarchy_counts()`

summarize hierarchy

#### Usage

    LFQDataSummariser$hierarchy_counts()

------------------------------------------------------------------------

### `LFQDataSummariser$hierarchy_counts_sample()`

number of elements at each level in every sample

#### Usage

    LFQDataSummariser$hierarchy_counts_sample(
      value = c("wide", "long"),
      nr_children = 1
    )

#### Arguments

- `value`:

  `"wide"` (one column per hierarchy level) or `"long"` format

- `nr_children`:

  get summary for 1,2 or more number of children

#### Returns

tibble of counts per sample and hierarchy level

------------------------------------------------------------------------

### `LFQDataSummariser$plot_hierarchy_counts_sample()`

barplot showing number of elements at each level in every sample

#### Usage

    LFQDataSummariser$plot_hierarchy_counts_sample(nr_children = 1)

#### Arguments

- `nr_children`:

  get summary for 1,2 or more number of children

------------------------------------------------------------------------

### `LFQDataSummariser$interaction_missing_stats()`

missing per condition and protein

#### Usage

    LFQDataSummariser$interaction_missing_stats()

------------------------------------------------------------------------

### `LFQDataSummariser$upset_interaction_missing_stats()`

upset plot with missing information per protein and condition

#### Usage

    LFQDataSummariser$upset_interaction_missing_stats(tr = 2)

#### Arguments

- `tr`:

  if less than tr observations in condition then missing

------------------------------------------------------------------------

### `LFQDataSummariser$missingness_per_group()`

missing stats per condition

#### Usage

    LFQDataSummariser$missingness_per_group()

------------------------------------------------------------------------

### `LFQDataSummariser$missingness_per_group_cumsum()`

missing stats per condition as cumulative sum

#### Usage

    LFQDataSummariser$missingness_per_group_cumsum()

------------------------------------------------------------------------

### `LFQDataSummariser$plot_missingness_per_group()`

barplot with number of features with 1,2, etc missing in condition

#### Usage

    LFQDataSummariser$plot_missingness_per_group()

#### Returns

ggplot

------------------------------------------------------------------------

### `LFQDataSummariser$plot_missingness_per_group_cumsum()`

barplot with cumulative sum of features with 1,2, etc missing in
condition

#### Usage

    LFQDataSummariser$plot_missingness_per_group_cumsum()

#### Returns

ggplot

------------------------------------------------------------------------

### `LFQDataSummariser$percentage_abundance()`

Does roll up to highest hierarchy and Computes the percent abundance of
proteins overall and within each group

#### Usage

    LFQDataSummariser$percentage_abundance()

#### Returns

data frame

------------------------------------------------------------------------

### `LFQDataSummariser$clone()`

The objects of this class are cloneable with this method.

#### Usage

    LFQDataSummariser$clone(deep = FALSE)

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

data <- istar$data
lfqdata <- LFQData$new(data, istar$config)
sum <- lfqdata$get_Summariser()
sum
#> <LFQDataSummariser>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     hierarchy_counts: function () 
#>     hierarchy_counts_sample: function (value = c("wide", "long"), nr_children = 1) 
#>     initialize: function (lfqdata) 
#>     interaction_missing_stats: function () 
#>     lfq: LFQData, R6
#>     missingness_per_group: function () 
#>     missingness_per_group_cumsum: function () 
#>     percentage_abundance: function () 
#>     plot_hierarchy_counts_sample: function (nr_children = 1) 
#>     plot_missingness_per_group: function () 
#>     plot_missingness_per_group_cumsum: function () 
#>     upset_interaction_missing_stats: function (tr = 2) 
sum$hierarchy_counts()
#> # A tibble: 1 × 2
#>   protein_Id peptide_Id
#>        <int>      <int>
#> 1         10         28
sum$hierarchy_counts_sample("wide")
#> # A tibble: 12 × 3
#>    sampleName protein_Id peptide_Id
#>    <chr>           <int>      <int>
#>  1 A_V1               10         23
#>  2 A_V2               10         26
#>  3 A_V3               10         25
#>  4 A_V4               10         23
#>  5 B_V1                8         26
#>  6 B_V2                9         26
#>  7 B_V3               10         27
#>  8 B_V4               10         26
#>  9 Ctrl_V1            10         24
#> 10 Ctrl_V2            10         25
#> 11 Ctrl_V3             9         24
#> 12 Ctrl_V4            10         25
sum$hierarchy_counts_sample("long")
#> # A tibble: 24 × 3
#>    sampleName key           nr
#>    <chr>      <chr>      <int>
#>  1 A_V1       protein_Id    10
#>  2 A_V1       peptide_Id    23
#>  3 A_V2       protein_Id    10
#>  4 A_V2       peptide_Id    26
#>  5 A_V3       protein_Id    10
#>  6 A_V3       peptide_Id    25
#>  7 A_V4       protein_Id    10
#>  8 A_V4       peptide_Id    23
#>  9 B_V1       protein_Id     8
#> 10 B_V1       peptide_Id    26
#> # ℹ 14 more rows
sum$plot_hierarchy_counts_sample()

sum$plot_hierarchy_counts_sample()

tmp <- sum$interaction_missing_stats()

sum$missingness_per_group()
#> # A tibble: 3 × 6
#>   group_ nrReplicates   `0`   `1`   `2`   `3`
#>   <chr>         <int> <int> <int> <int> <int>
#> 1 A                 4    19     5     2     2
#> 2 B                 4    21     7    NA    NA
#> 3 Ctrl              4    17     8     3    NA
sum$missingness_per_group_cumsum()
#> # A tibble: 3 × 6
#>   group_ nrReplicates   `0`   `1`   `2`   `3`
#>   <chr>         <int> <int> <int> <int> <int>
#> 1 A                 4    19    24    26    28
#> 2 B                 4    21    28    NA    NA
#> 3 Ctrl              4    17    25    28    NA
sum$plot_missingness_per_group()

sum$plot_missingness_per_group_cumsum()

sum$upset_interaction_missing_stats()

sum$percentage_abundance()
#> # A tibble: 112 × 16
#>    protein_Id  peptide_Id group_ nrReplicates nrMeasured nrNAs    sd   var
#>    <chr>       <chr>      <chr>         <int>      <int> <int> <dbl> <dbl>
#>  1 BEJI92~5282 HBkZvdhT   All              12          8     4  1.35  1.82
#>  2 DoWup2~5896 KVUnZ6oZ   All              12         12     0  2.55  6.48
#>  3 Fl4JiV~8625 KpyeEoiy   All              12         10     2  1.30  1.68
#>  4 0EfVhX~0087 dJkdz7so   All              12          9     3  4.36 19.0 
#>  5 Fl4JiV~8625 fv2Ck8hz   All              12          8     4  4.31 18.6 
#>  6 9VUkAq~4703 eIC06D7g   All              12         11     1  4.97 24.7 
#>  7 Fl4JiV~8625 GsUIOl6Q   All              12         10     2  3.29 10.8 
#>  8 HvIpHG~9079 opjydeWJ   All              12         10     2  3.16  9.99
#>  9 SGIVBl~5782 03kMlNn1   All              12          9     3  1.24  1.53
#> 10 BEJI92~5282 qQ1GK8Un   All              12         11     1  3.69 13.6 
#> # ℹ 102 more rows
#> # ℹ 8 more variables: meanAbundance <dbl>, medianAbundance <dbl>, CV <dbl>,
#> #   interaction <chr>, id <int>, abundance_percent <dbl>,
#> #   abundance_percent_cumulative <dbl>, percent_prot <dbl>
```
