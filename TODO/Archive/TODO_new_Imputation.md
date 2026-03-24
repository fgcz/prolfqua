# New Imputation

I want an alternative lm based implementaton of ContrastImpute, I would be an extension of the contrast lm class.                                                                          

We build the models, as usual. Then for those proteins where we fail to build, we impute the NA with the lod. but also set all values below lod to lod. We then fit the model, however we do not use the convariance matrix, of the model fit. But we use a covariance matrix, which is a mean or median of the covariance matrices which we did estimate from the other data.

The reason being, since we impute with a constant value we would underestimate the variance, therefore we take the variance from the other models.

The idea of using the lod not just replace the NA but also as minimum value is to avoid situations that a group where we observed the data will get a higher estimate that a group where we failed to obtain an estimate.

Last but not least degrees of freedom, ideally only as many as indeed observed values.

So basic principle impute -> fit model -> update model to reflect missingness and observed variances. What do you think (no coding here. we discuss)



Interesting approach. Let me think through the mechanics and flag some considerations.

### What I like

**Reusing the lm infrastructure** — you get contrast matrices, multi-factor designs, and the existing `Contrasts` class machinery for free. No manual group-mean arithmetic.

**Clamping to LOD (not just imputing NAs)** — this solves the current zero-out logic problem elegantly. If group A has one observation at 8.5 and group B is all-missing (imputed at LOD=10), the current code gives a "wrong direction" fold change and zeroes it out. With clamping, group A's value becomes max(8.5, LOD) and the contrast naturally behaves.

**Borrowing the covariance matrix** — this is the key insight. The imputed fit gives you correct *point estimates* (coefficients), but the variance is artificial. Swapping in a typical covariance matrix from successful fits restores honest uncertainty.

### Questions and concerns

**1. Which covariance matrix to borrow?**

You said mean or median of covariance matrices from successful fits. But covariance matrices aren't scalars — you need to aggregate them sensibly:

- **Option A:** Take the median of `sigma^2` (residual variance) across successful fits, then scale the `(X'X)^{-1}` from the imputed fit by that median sigma. This is simple and preserves the design structure.
- **Option B:** Take the full `vcov()` matrices and compute element-wise median. This is messier and the matrices may not be positive definite after aggregation.



Option A seems cleaner — you're essentially saying "the coefficients are from the imputed fit, but the residual variance is borrowed." This is conceptually close to what limma's `squeezeVar` does.

COMMENT: let's have both options A and B, build_model has a parameter to decide which.

**2. Degrees of freedom — which df?**

You want df to reflect actual observed values. Two choices:

- **Residual df of the imputed model** (n_total - p) — but this overcounts since imputed values aren't real observations
- **n_observed - p** — honest but could be very small or even zero (which is why the model failed in the first place)

If n_observed < p, you have no real df. You could use the df from the borrowed variance estimate (the successful fits), which gives you a "prior df" similar to limma's approach.

COMMENT: again, two options max(n_observed - p, 1) or df from the borrowed, again param.

**3. When exactly do you impute?**

You said: fit normally first, then impute for failures. But the `strategy_lm()` / `build_model()` path fits per-protein. So:

- Proteins with enough data → normal fit, normal covariance → use as-is (these are "donors")
- Proteins with too many NAs → model fails → impute with LOD + clamp → refit → replace covariance with donor estimate

This creates two tiers of results. The interface should probably flag which tier each protein came from (similar to the current `indic` column).
Where do we have the indic column 

COMMENT: we have the model column, one is e.g. WaldTEst the other WaldTestImputed or something. right? Same way as merge_contrast_results... howewer, I guess, we do not need to fit for all but only those with missing.



**4. What does "model fails" mean exactly?**

- `lm()` returns NA coefficients? 
- Singular fit?
- Or a threshold like "more than X% missing"?

COMMENT: For the start, we do not focus on model, we focus on the contrasts, if not estimated we must refit the protein, and reistmate the contrast, or we can focus on the models, basically if any of the coefficients is NA, we will fails in estimating contrasts. So answer yes NA coeffs , if NA coeffs then singular? or can it be that NA coefs and not singular? So basically we would have Model class, with a modelDF with some additional annotations. 


You might also want to impute for proteins that *technically* fit but have, say, only 1 observation in one group — the fit succeeds but the variance estimate is garbage. COMMENT: This is a good point! I believe limma requires 2 obs per group, while our lm fit only one. 




**5. Relationship to existing classes**

Architecturally, this could be:

- A **post-processor on `Model`** — fit all proteins (imputing where needed), then selectively replace `vcov` on the model objects.


- A **new `Contrasts*` class** that wraps `Contrasts` — takes a normal `Model`, identifies failures, refits with imputation, borrows covariance, and merges results

COMMENT: would prefer, that build_model handles all this, and Contrasts stays unchanged.

- A **facade** like `ContrastsRLMFacade` that orchestrates: fit → identify failures → impute → refit → patch covariance → merge

COMMENT: We can have a separate facade, but all it does is to parameterize build_models with impute.

The facade pattern seems most natural given how the other "fallback" approaches work (e.g., `ContrastsMissing` is already used as a fallback in `prolfquapp`). COMMENT: disagree in this case. ContrastMissing is an exception, since we have to build two contrasts and merge.




### Summary of the design as I understand it

```
For each protein:
  1. Try lm fit on observed data
  2. If success → keep fit, contribute vcov to "donor pool"
  3. If failure:
     a. Replace NAs with LOD, clamp all values to max(value, LOD)
     b. Fit lm on imputed data → get coefficients
     c. Replace vcov with donor estimate (median sigma² × (X'X)⁻¹)
     d. Set df = n_observed - p (or donor df if too small)
  4. Extract contrasts using standard Contrasts machinery
```

Does this match your thinking? The main open question from my side is whether option A (borrow scalar sigma², keep imputed design's `(X'X)^{-1}`) is what you have in mind for the covariance replacement.