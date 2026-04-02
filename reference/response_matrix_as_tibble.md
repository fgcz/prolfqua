# Takes matrix of responses and converts into tibble

Takes matrix of responses and converts into tibble

## Usage

``` r
response_matrix_as_tibble(pdata, value, config, data = NULL, sep = "~lfq~")
```

## Arguments

- pdata:

  (matrix)

- value:

  name of column to store values in. (see \`gather\`)

- config:

  AnalysisConfiguration

- data:

  lfqdata

- sep:

  separater to unite the hierarchy keys.

## Examples

``` r
dd <- prolfqua::sim_lfq_data_peptide_config()
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
data <- dd$data
conf <- dd$config
res <- tidy_to_wide_config(data, conf, as.matrix = TRUE)

res <- scale(res$data)
xx <- response_matrix_as_tibble(res,"srm_intensityScaled", conf)
xx <- response_matrix_as_tibble(res,"srm_intensityScaled", conf, data)
#> Joining with `by = join_by(sampleName, isotopeLabel, protein_Id, peptide_Id)`
conf$get_response() == "srm_intensityScaled"
#> [1] TRUE
```
