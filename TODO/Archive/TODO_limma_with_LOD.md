we did implement lm with imputation see @DONE_LM_Impute.md 
based on build_model_impute and the default Contrast R 6 class.

I would like to extend the limma in the same fashion.

Two options:

Option 1:
a) use data with NA fit limma,
b) use data with LOD fit limma, 
determine failed models in a), and replace those with the fits from b).

or just run a) and augment failed models similar to our lm_missing build.


Now I have a clear picture. Here's my analysis of the TODO:

---

The TODO asks about extending limma with LOD imputation, analogous to what `build_model_impute` does for `lm`. There are two approaches mentioned:

**Option 1 (fit twice):** Fit limma on NA data, fit limma on LOD-imputed data, replace failed proteins from (a) with fits from (b).

**Option 2 (augment failed, like `lm_missing`):** Fit limma on NA data, then augment proteins that failed (all-missing in one group) with `ContrastsMissing`.

The key difference from `lm` is that **limma is a matrix-based pipeline** — `lmFit` fits all proteins at once, not per-protein. This changes the mechanics:

**For Option 1 (limma + LOD imputation):**
- `lmFit` on NA data → proteins with entire groups missing get `NA` coefficients
- `lmFit` on LOD-imputed data → all proteins get coefficients, but variance is underestimated for imputed ones
- Challenge: you can't easily "borrow sigma" per-protein like `build_model_impute` does, because limma's `eBayes` pools variance across all proteins. You'd need to either:
  - Replace individual rows in `fit$coefficients` / `fit$sigma` from the imputed fit (messy, breaks limma's internal assumptions)
  - Or keep the two fits separate and merge at the contrast result level (cleaner)

**For Option 2 (limma + ContrastsMissing merge):**
COMMENT: most stupid idea I ever read. Contrast Missing has sever limitation.

- This is simpler and mirrors `ContrastsLMMissingFacade` exactly
- Fit limma on NA data → `ContrastsLimma` → merge with `ContrastsMissing`
- Proteins with enough data get limma results; completely-missing-group proteins get fold-change-from-missingness estimates
- This is essentially already possible with the existing pieces

COMMENT: you guid must be CONTRAST_LM_IMPUTE! look how we transfer variance and degrees of freedom! It is not about easy it is about improving.

**My assessment:**
- **Option 1 is more ambitious** but the merge-at-result-level variant is feasible: run `ContrastsLimma` on NA data, run `ContrastsLimma` on LOD-imputed data (with `eBayes = FALSE` + borrowed sigma), and merge results — replacing failed proteins from the first with results from the second.

Comment: explain  (with `eBayes = FALSE` + borrowed sigma).

