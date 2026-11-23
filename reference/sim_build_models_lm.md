# Build linear models on simulated data for testing

Build linear models on simulated data for testing

## Usage

``` r
sim_build_models_lm(
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

a [Model](https://fgcz.github.io/prolfqua/reference/Model.md) object
with one linear model per simulated protein, fitted with
[`build_model()`](https://fgcz.github.io/prolfqua/reference/build_model.md)

## Note

Deprecated in prolfqua 1.9.5; use
`sim_lfq_data_protein_config() with build_contrast_analysis()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
