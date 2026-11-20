# Compute fold changes given Contrasts

Compute fold changes given Contrasts

## Usage

``` r
get_contrast(data, hierarchy_keys, contrasts)
```

## Arguments

- data:

  hierarchy_keys of Analysis Configuration

- contrasts:

  list of contrasts to compute

## See also

Other imputation:
[`missigness_histogram()`](https://fgcz.github.io/prolfqua/reference/missigness_histogram.md),
[`missingness_per_condition()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://fgcz.github.io/prolfqua/reference/missingness_per_condition_cumsum.md),
[`upset_interaction_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_interaction_missing_stats.md),
[`upset_missing_stats()`](https://fgcz.github.io/prolfqua/reference/upset_missing_stats.md)

## Examples

``` r

istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfq <- LFQData$new(istar$data, istar$config)
lfq$complete_cases()

Contrasts <- c("dilution.b-a" = "group_A - group_B", "dilution.c-e" = "group_A - group_Ctrl")

var <- summarize_stats(lfq)
var <- prolfqua::make_interaction_column(var, columns = lfq$relevant_factor_keys())

imp <- var |> tidyr::pivot_wider(id_cols = lfq$hierarchy_keys(),
                        names_from = interaction,
                        values_from = !!rlang::sym("meanAbundance"))

imputed <- get_contrast(imp, lfq$hierarchy_keys(), Contrasts)

```
