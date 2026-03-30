# Convert LFQData to limma EList with design matrix and metadata

Shared preamble for all `build_model_limma*` and `build_model_limpa`
functions. Pivots LFQData to wide format, builds the design matrix from
the formula, resolves the subject_Id / isotopeLabel, and creates a dummy
lm for linfct extraction.

## Usage

``` r
.lfqdata_to_elist(lfqdata, formula)
```

## Arguments

- lfqdata:

  an [`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md)
  object

- formula:

  a formula (with response and RHS)

## Value

a list with components:

- elist:

  limma EList with `$E` = expression matrix

- expr_matrix:

  the expression matrix (same as `elist$E`)

- design:

  the design matrix

- annotation:

  sample-level annotation data.frame

- subject_Id:

  character vector of hierarchy keys (possibly including isotopeLabel)

- rowdata:

  data.frame with one row per feature, columns = subject_Id

- rhs_formula:

  the RHS-only formula

- dummy_model:

  a dummy `lm` fitted on one complete row
