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
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
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
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id, peptide_Id)`
```
