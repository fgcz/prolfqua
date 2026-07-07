# Vectorized contrast path diverges from loop path on rank-deficient models

**Date:** 2026-07-06
**Status:** FIXED (2026-07-06). Guard at `R/tidyMS_contrasts.R:541` now counts nonzero weights;
deterministic regression test added at `tests/testthat/test-vectorize-contrasts.R`
("rank-deficient model, canceling NA weights"); NEWS.md updated. Loop and vectorized paths verified
to both return `NA` for the reproduction below.
**Severity:** Moderate — only triggers under `options(prolfqua.vectorize = TRUE)` (non-default) *and* a
rank-deficient (aliased) model fit. Default runs (`prolfqua.vectorize = FALSE`) are unaffected.

## Problem

`compute_contrast_vectorized()` (`R/tidyMS_contrasts.R:541`) and the loop-based `compute_contrast()`
are documented as numerically identical, but they diverge for a contrast row whose weights fall on
missing (`NA`, i.e. aliased / rank-deficient) coefficients:

- **Loop path** (`compute_contrast`): correctly returns `estimate = NA, p.value = NA` — it refuses to
  evaluate a contrast that references coefficients not estimable in the model.
- **Vectorized path**: returns `estimate = 0, p.value = NaN` — it treats the row as valid and the
  zero-filled coefficients produce a spurious fold-change of 0.

A fold-change of exactly 0 with a `NaN` p-value can be misread downstream as a genuine "no change"
result rather than an undefined contrast.

## Root cause

Two steps interact. First, the estimate is **not** computed as `1*NA + (-1)*NA` (which would be
`NA`) — the vectorized code deliberately zero-fills the `NA` coefficients before the matrix multiply
(`R/tidyMS_contrasts.R:535-537`):

```r
coef_zero <- coef_aligned
coef_zero[na_coefs] <- 0            # NA coefficients replaced by 0 BEFORE the multiply
estimate <- as.vector(linfct %*% coef_zero)
```

So for a `+1 / -1` contrast on two `NA` coefficients the raw estimate is `1*0 + (-1)*0 = 0`, not
`NA`. (The zero-fill exists so that columns with genuinely zero weight don't poison the estimate,
since `0 * NA` is `NA` in R.)

The estimate is therefore always numeric; the *only* thing that restores `NA` for rank-deficient
rows is the guard on `R/tidyMS_contrasts.R:541`:

```r
invalid <- as.logical(linfct[, na_coefs, drop = FALSE] %*% rep(1, sum(na_coefs)) != 0)
```

This flags a row invalid only when the **signed sum** of its weights on the `NA` columns is nonzero.
For `+1` and `-1` the signed sum is `1 + (-1) = 0`, so `invalid = FALSE`, the guard does not fire,
and the spurious `0` estimate is kept. Verified intermediate values: `coef_zero` on the NA columns =
`0 0`; raw `estimate = 0`; signed-sum guard = `0`; `invalid = FALSE`. The fix below (counting
nonzero weights: `|1| + |-1| = 2`) makes `invalid = TRUE`.

## Reproduction

```r
suppressMessages(devtools::load_all("."))
set.seed(1)
d <- data.frame(
  y = rnorm(12),
  A = factor(rep(c("a1","a2","a3"), each = 4)),
  B = factor(rep(c("b1","b2","b3"), each = 4))   # B perfectly aliased with A
)
m <- lm(y ~ A + B, data = d)                      # Bb2, Bb3 coefficients are NA (aliased)
na_names <- names(which(is.na(coefficients(m))))
cn <- names(coefficients(m))
L <- matrix(0, nrow = 1, ncol = length(cn), dimnames = list("cancel", cn))
L[1, na_names[1]] <- 1
L[1, na_names[2]] <- -1                            # +1/-1 on two NA coefficients -> signed sum 0

options(prolfqua.vectorize = FALSE); o <- compute_contrast(m, L)
options(prolfqua.vectorize = TRUE);  v <- compute_contrast_vectorized(m, L)
options(prolfqua.vectorize = FALSE)
# o$estimate == NA, o$p.value == NA   (correct)
# v$estimate == 0,  v$p.value == NaN  (bug)
```

## Fix

Detect *any* nonzero weight on an `NA` coefficient instead of summing signed weights
(`R/tidyMS_contrasts.R:541`):

```r
invalid <- as.logical((linfct[, na_coefs, drop = FALSE] != 0) %*% rep(1, sum(na_coefs)) != 0)
```

## Why this wasn't caught

The equivalence test for NA-coefficient models (`tests/testthat/test-vectorize-contrasts.R:245`,
which asserts `is.na(orig$estimate) == is.na(vec$estimate)`) is **skipped** at runtime:
`skip_if(!any(has_na), "No models with NA coefficients in simulated data")` — the simulated data
never produces a rank-deficient fit, so the check never runs.

When fixing, replace the skipped, data-dependent path with a hand-built aliased `lm()` (as in the
reproduction above) so the loop-vs-vectorized equivalence is actually exercised for rank-deficient
models. Add a `NEWS.md` entry under the current version.
