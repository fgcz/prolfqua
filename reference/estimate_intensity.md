# Aggregates e.g. protein abundances from peptide abundances

Aggregates e.g. protein abundances from peptide abundances

## Usage

``` r
estimate_intensity(lfqdata, method = c("medpolish", "rlm"), control = list())
```

## Arguments

- lfqdata:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- method:

  "medpolish" (Tukey's median polish, column `medpolish`) or "rlm"
  (robust regression with
  [`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html), column `lmrob`)

- control:

  convergence controls merged over the defaults. medpolish: `maxiter`
  (default 10) and `eps` (default 0.01; the polish stops when the sum of
  absolute residuals changes by less than `eps` times that sum, a
  relative criterion independent of the input scale); rlm: `maxit`
  (default 20) and `acc` (default 1e-4, relative change of the IRLS
  residuals).

## Value

returns list with data (data.frame), config (AnalysisConfiguration) and
failures (tibble with the subject `id` and `warning` text of every
warning raised while estimating one subject; the warnings are not
re-raised)

## Note

Deprecated in prolfqua 1.9.5; use
`LFQData$get_Aggregator("medpolish")$aggregate()`. See
[prolfqua-deprecated](https://fgcz.github.io/prolfqua/reference/prolfqua-deprecated.md).
