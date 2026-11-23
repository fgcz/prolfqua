# Fit a mixed-effects model to simulated data for examples

Fit a mixed-effects model to simulated data for examples

## Usage

``` r
sim_make_model_lmer(
  model = c("parallel2", "parallel3", "factors", "interaction")
)
```

## Arguments

- model:

  design: `"parallel2"`, `"parallel3"`, `"factors"` or `"interaction"`

## Value

an `lmerModLmerTest` object, the first successful fit among ten
simulated proteins

## Note

Deprecated in prolfqua 1.9.5; use
`sim_lfq_data_peptide_config() with build_contrast_analysis()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
