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
[`UpSet_interaction_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_interaction_missing_stats.md),
[`UpSet_missing_stats()`](https://wolski.github.io/prolfqua/reference/UpSet_missing_stats.md),
[`missigness_histogram()`](https://wolski.github.io/prolfqua/reference/missigness_histogram.md),
[`missingness_per_condition()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition.md),
[`missingness_per_condition_cumsum()`](https://wolski.github.io/prolfqua/reference/missingness_per_condition_cumsum.md)

## Examples

``` r

istar <- sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
config <- istar$config
analysis <- istar$data
data <- complete_cases(analysis, config)
#> completing cases

Contrasts <- c("dilution.b-a" = "group_A - group_B", "dilution.c-e" = "group_A - group_Ctrl")

var = summarize_stats(data, config)
#> completing cases
var <- prolfqua::make_interaction_column(var, columns = config$factor_keys_depth())

imp <- var |> tidyr::pivot_wider(id_cols = config$hierarchy_keys(),
                        names_from = interaction,
                        values_from = !!rlang::sym("meanAbundance"))

imputed <- get_contrast(imp, config$hierarchy_keys(), Contrasts)
#> dilution.b-a=group_A - group_B
#> dilution.c-e=group_A - group_Ctrl

```
