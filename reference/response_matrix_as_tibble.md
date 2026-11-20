# Takes matrix of responses and converts into tibble

Takes matrix of responses and converts into tibble

## Usage

``` r
response_matrix_as_tibble(pdata, value, config, data = NULL, sep = "~lfq~")
```

## Arguments

- pdata:

  matrix with rownames encoding hierarchy keys

- value:

  name of column to store values in

- config:

  AnalysisConfiguration (needed for column name mapping)

- data:

  optional data.frame to join back to

- sep:

  separator used to unite the hierarchy keys

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- prolfqua::LFQData$new(dd$data, dd$config)
res <- tidy_to_wide_config(lfqdata, as.matrix = TRUE)

res_scaled <- scale(res$data)
xx <- response_matrix_as_tibble(
  res_scaled, "srm_intensityScaled", lfqdata$get_config()
)
xx <- response_matrix_as_tibble(
  res_scaled, "srm_intensityScaled",
  lfqdata$get_config(), lfqdata$data_long()
)
```
