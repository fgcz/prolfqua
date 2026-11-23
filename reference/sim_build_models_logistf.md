# Build Firth logistic models on simulated data for testing

Build Firth logistic models on simulated data for testing

## Usage

``` r
sim_build_models_logistf(
  model = c("parallel2", "parallel3", "factors", "interaction"),
  n_proteins = 10,
  with_missing = TRUE,
  weight_missing = 1,
  Nprot
)
```

## Arguments

- model:

  design: `"parallel2"`, `"parallel3"`, `"factors"` or `"interaction"`

- n_proteins:

  number of simulated proteins

- with_missing:

  include simulated missing measurements

- weight_missing:

  weight controlling simulated missingness

- Nprot:

  deprecated; use `n_proteins`

## Value

a [ModelFirth](https://fgcz.github.io/prolfqua/reference/ModelFirth.md)
object fitted with
[`build_model_logistf()`](https://fgcz.github.io/prolfqua/reference/build_model_logistf.md)

## Note

Deprecated in prolfqua 1.9.5; use
`sim_lfq_data_protein_config() with build_contrast_analysis()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
