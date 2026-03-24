# How `lmFit` Works Internally

`lmFit` is essentially a **dispatcher**. It extracts the expression matrix via `getEAWP()`, validates the design matrix, and then routes to one of three lower-level fitting functions depending on the arguments:

1. **`lm.series`** — ordinary (weighted) least squares. Used when `method="ls"` and there's no blocking/duplicate correlation structure. This is the most common path.
2. **`gls.series`** — generalized least squares. Used when there are duplicate probes (`ndups > 1`) or a `block` factor with a known `correlation`.
3. **`mrlm`** — robust regression via `MASS::rlm`. Used when `method="robust"`.

All three return an `MArrayLM` object.

---

## The Core Fitting in `lm.series` (the Default Path)

This is where the magic is. There are **two branches**:

### Fast Path (no per-gene missing values or per-probe weights)

The design matrix QR decomposition is the **same for every gene**, so limma does it *once* and fits all genes simultaneously via a single call to R's `lm.fit(design, t(M))`. This is why limma is so fast — it's a single LAPACK/BLAS call for thousands of genes at once, not a loop.

The covariance factor is computed as:

```r
fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
```

This is `(X'X)^{-1}` (or the weighted equivalent), computed once from the Cholesky/QR factorization of the design. It is **shared across all genes** because the design is the same.

The per-gene residual standard error `sigma[g]` is estimated from:

```r
fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank + 1):narrays, , drop=FALSE]^2))
```

### Slow Path (missing values or per-probe weights)

When genes have different patterns of missing data or individual observation weights, each gene needs its own QR decomposition, so limma falls back to a `for` loop calling `lm.fit` or `lm.wfit` gene by gene. At the end, a single `cov.coefficients` is still computed from the *full* design matrix.

---

## The `gls.series` Path

For correlated observations (duplicate spots or blocked samples), limma constructs a correlation matrix `V`, takes its Cholesky decomposition, and "decorrelates" both `Y` and `X` by back-solving:

```r
y <- backsolve(cholV, t(M), transpose = TRUE)
X <- backsolve(cholV, design, transpose = TRUE)
fit <- lm.fit(X, y)
```

This transforms the GLS problem into an OLS problem on the whitened data, then proceeds identically to `lm.series`.

---

## Accessing the Covariance Matrix

After fitting:

```r
fit <- lmFit(y, design)
```

The object `fit` contains:

- **`fit$cov.coefficients`** — this is `(X'WX)^{-1}`, the **unscaled** covariance matrix of the coefficients. It's a single `p × p` matrix (where `p` = number of design columns) that is **shared across all genes**.
- **`fit$sigma`** — a vector of per-gene residual standard deviations (one per gene/protein).
- **`fit$stdev.unscaled`** — a matrix (genes × coefficients) containing `sqrt(diag(cov.coefficients))`, i.e., the unscaled standard errors per coefficient.

### Reconstructing the Full Per-Gene Covariance Matrix

The full covariance matrix of the estimated coefficients for gene `g` is:

```r
# Per-gene covariance of beta_g:
cov_g <- fit$sigma[g]^2 * fit$cov.coefficients
```

This works because limma assumes the same design matrix for all genes — the only thing that differs per gene is `sigma^2`. So the covariance structure `(X'X)^{-1}` is shared, and you just scale it by the gene-specific variance estimate.

After `eBayes()`, you'd use the **moderated** (shrunken) variance instead:

```r
fit <- eBayes(fit)
cov_g_moderated <- fit$s2.post[g] * fit$cov.coefficients
```

where `fit$s2.post` is the posterior variance after the empirical Bayes shrinkage.

### Important Caveat

If some genes have missing values (so they were fit on different subsets of samples), the `cov.coefficients` matrix stored in the fit object is computed from the **full design matrix** and may not be exactly right for genes with missing data. For those genes, the per-gene `stdev.unscaled` values (which *are* computed gene-by-gene in the slow path) are more reliable than reconstructing from the shared `cov.coefficients`. In practice, this matters only when you have substantial missingness.

---

## Quick Reference Table

| What you want                | How to get it                              |
|------------------------------|--------------------------------------------|
| Unscaled covariance (shared) | `fit$cov.coefficients`                     |
| Per-gene σ²                  | `fit$sigma^2`                              |
| Full per-gene Cov(β̂)        | `fit$sigma[g]^2 * fit$cov.coefficients`    |
| After eBayes moderation      | `fit$s2.post[g] * fit$cov.coefficients`    |
| Per-gene SE of each coef     | `fit$stdev.unscaled * fit$sigma`           |

The covariance information is there — it's just stored in a factored form (shared structure × per-gene scalar) rather than as thousands of separate matrices, which is both memory-efficient and mathematically natural for the parallel regression framework.

---

# Deep Dive: `gls.series`, Duplicate Probes, Blocking, and Paired Designs

## Duplicate Probes (`ndups`) vs. Blocking (`block`)

These are two **different** correlation structures that `gls.series` can handle:

### Duplicate Probes (`ndups`)

This is a microarray-specific concept. On older spotted arrays, the same gene/probe was physically printed multiple times on a single slide — e.g., every gene appears twice. The `ndups=2, spacing=1` parameters describe the physical layout (consecutive rows = same probe). These within-array replicates are correlated because they come from the same hybridization on the same sample. This is *not* about biological replication at all — it's about technical replication within a single chip.

### Blocking (`block`)

This is the more general and more relevant concept for modern experiments. This is where you have **groups of correlated samples across arrays**. A paired experiment is the classic example: tumor and normal tissue from the same patient. The two arrays from patient 1 are more correlated with each other than with arrays from patient 2, because they share patient-level biology.

**`block` is how you handle paired designs in limma**, but it's more general than just pairs — it covers any grouping that induces correlation (e.g., repeated measures, samples from the same batch, longitudinal timepoints from the same subject).

---

## Is It the Same as Averaging Over Duplicates?

No, and this matters. Averaging duplicate probes (which you can do with `avereps()` in limma) is a simpler approach, but it:

- Discards information about within-duplicate variability
- Gives you correct point estimates (the coefficients are the same)
- But gives you **wrong standard errors** — they'll typically be too small because you've artificially reduced the apparent number of observations without accounting for the fact that they weren't independent

The GLS approach properly models the correlation. It says "these observations are not independent, and I'll account for that when computing standard errors and test statistics." The Cholesky whitening trick in the source code (`backsolve(cholV, ...)`) is exactly this — it transforms the correlated data into uncorrelated data with correct variance structure, then fits OLS on the transformed data.

That said, for duplicate *probes* specifically, `avereps()` followed by `lm.series` is a reasonable practical shortcut that many people use, and the results are usually very similar. For paired/blocked *samples*, you really should use the proper GLS approach.

---

## How Is This Encoded?

This is a key point: **the correlation structure is NOT in the design matrix**. Limma separates the two concerns:

### Design Matrix → Fixed Effects

The design matrix handles what groups exist and what comparisons you want. For a paired tumor/normal experiment:

```r
# The design matrix encodes the treatment effect
design <- model.matrix(~ condition)
# e.g., columns: Intercept, conditionTumor
```

### Block + Correlation → Dependence Structure

These are specified as *separate arguments*:

```r
# block encodes which observations are correlated
block <- patient_id   # e.g., c(1,1,2,2,3,3,...)

# First estimate the consensus correlation
corfit <- duplicateCorrelation(y, design, block = patient_id)

# Then pass both to lmFit
fit <- lmFit(y, design, block = patient_id,
             correlation = corfit$consensus.correlation)
```

The `correlation` parameter is a **single scalar** — limma assumes a compound symmetry structure (equal correlation between all observations in the same block). This is estimated across all genes by `duplicateCorrelation()`, which fits a mixed model with a random intercept per block and then takes a consensus (trimmed mean) of the per-gene correlation estimates.

---

## The Alternative: Encoding Pairing in the Design Matrix

You *could* handle a paired design purely through the design matrix by including patient as a fixed effect:

```r
design <- model.matrix(~ patient + condition)
```

This is perfectly valid statistically and doesn't require `block`/`correlation` at all — `lm.series` handles it directly. The tradeoff is:

- **Design matrix approach**: each patient gets its own coefficient, eating up degrees of freedom. With 50 patients you lose 49 df. But each gene gets its own patient effects.
- **Block/correlation approach**: estimates a single shared correlation parameter across all genes (empirical Bayes-like borrowing of information), preserves degrees of freedom, but assumes compound symmetry and a common correlation across genes.

For small experiments (few patients), the design matrix approach can cost you too many degrees of freedom, and the block approach is better. For large experiments, the design matrix approach is fine and avoids the compound symmetry assumption.

---

## Summary: Correlation Structures in limma

| Concept | What it handles | How it's specified | Statistical mechanism |
|---|---|---|---|
| `ndups`/`spacing` | Same probe printed multiple times on one array | Arguments to `lmFit` | GLS via Cholesky whitening |
| `block` | Correlated samples (paired, repeated measures) | Argument to `lmFit` + `duplicateCorrelation` | GLS via Cholesky whitening |
| Patient in design matrix | Paired design as fixed effects | Columns in design matrix | OLS with extra covariates |

---

# Concrete `ndups` Example and the Proteomics Peptide Problem

## What `ndups` Actually Looks Like

On an old spotted cDNA array, a gene might be physically printed twice in consecutive rows. The raw data matrix would look like this:

```
             Sample1  Sample2  Sample3
GeneA_spot1    3.2      4.1      3.8     ← row 1
GeneA_spot2    3.4      3.9      3.7     ← row 2
GeneB_spot1    5.1      5.3      4.9     ← row 3
GeneB_spot2    4.8      5.0      5.1     ← row 4
```

Every gene has **exactly 2** rows, in a perfectly regular pattern. You'd call:

```r
fit <- lmFit(y, design, ndups = 2, spacing = 1,
             correlation = corfit$consensus.correlation)
```

Internally, `unwrapdups()` reshapes this into a wider matrix where each gene has one row but doubled columns (the two spots become two "pseudo-arrays"), and the correlation between spots from the same gene on the same array is modeled via GLS.

The critical constraint: **`ndups` must be a single fixed integer, the same for every gene**, and the layout must be perfectly regular. This is fine for a manufactured array where every gene was printed exactly twice. It is completely unsuitable for variable numbers of peptides per protein.

---

## Why This Breaks Down for Peptides per Protein

In proteomics, you have something like:

```
             Sample1  Sample2  Sample3
ProteinA_pep1   3.2     4.1      3.8
ProteinA_pep2   3.4     3.9      3.7
ProteinA_pep3   2.9     3.5      3.3
ProteinB_pep1   5.1     5.3      4.9
ProteinB_pep2   4.8     5.0      5.1
ProteinC_pep1   6.2     6.0      6.4
ProteinC_pep2   6.5     6.3      6.1
ProteinC_pep3   6.0     5.8      6.2
ProteinC_pep4   6.3     6.1      5.9
```

Protein A has 3 peptides, B has 2, C has 4. There's no single `ndups` value that works here — and `unwrapdups()` would produce garbage.

---

## Actual Options for Variable Peptides per Protein

There are fundamentally two philosophies:

### 1. Summarize First, Then Test (Most Common)

Aggregate peptides → protein-level intensities, then run standard `lmFit` on the protein summaries. The aggregation method matters:

- **Median polish** (Tukey's method, used in MSstats) — iteratively subtracts row and column medians, robust to outlier peptides
- **maxLFQ** (used in MaxQuant, DIA-NN) — finds the protein intensity that best explains all pairwise peptide ratios between samples
- **Simple mean/median** of peptide log-intensities — quick and dirty, works surprisingly well in practice

This is what most proteomics pipelines do. The downside is that you **lose information about peptide-level variability** — a protein with 50 peptides gets the same weight as one with 2 peptides, even though you're much more certain about the first one.

You *can* partially recover this by passing `weights` to `lmFit`, e.g., weighting each protein by its number of peptides or by the inverse variance of its peptide-level estimates.

### 2. Model at Peptide Level with Mixed Effects (Statistically Better, More Complex)

This is what **msqrob2** does. It fits per-protein mixed models on the peptide-level data:

```
log_intensity_ij = protein_effect_i + peptide_j + treatment + error
```

where peptide is a random effect nested within protein. Each protein has its own model that properly accounts for however many peptides it has. Then it borrows the empirical Bayes shrinkage idea from limma for the variance estimation across proteins.

In R this looks like:

```r
# msqrob2 approach (Bioconductor)
library(msqrob2)
pe <- msqrob(pe, formula = ~ condition, modelColumnName = "msqrobModels")
```

Under the hood, for each protein it fits something like an `lmer()` model with peptide as a random intercept, then does the eBayes-style moderation across proteins.

### 3. The `block` Hack (Possible but Not Ideal)

You *could* try to use limma's `block` argument with protein as the blocking factor, treating each peptide as an "observation":

```r
corfit <- duplicateCorrelation(peptide_matrix, design, block = protein_id)
fit <- lmFit(peptide_matrix, design, block = protein_id,
             correlation = corfit$consensus.correlation)
```

This technically runs, but it has real problems. It estimates a **single consensus correlation** across all proteins, which makes little sense when proteins range from 1 to 1000 peptides. And it treats each peptide as having equal information about the treatment effect, ignoring that some peptides are noisy or have missing values in different patterns. It also conflates the peptide-level model (which should have a peptide-specific intercept) with the sample-level correlation structure.

---

## Practical Recommendations for Proteomics

| Approach | Pros | Cons |
|---|---|---|
| **Summarize (maxLFQ / median polish) → limma** | Simple, fast, well-understood, battle-tested | Loses peptide-level variance info |
| **Summarize → limma with peptide-count weights** | Simple improvement over above | Weights are approximate |
| **msqrob2** | Statistically principled, handles variable peptides, proper mixed models + eBayes | More complex, slower, less widely used |
| **limma `block` on peptide data** | Uses existing limma machinery | Single correlation is a poor assumption here |

Most published proteomics studies use option 1 and it works well. If you care about properly propagating uncertainty from the peptide level (especially with lots of missingness), msqrob2 is the right tool — it was designed precisely for this problem.

---

# Why Constant Per-Protein Weights Do Nothing in limma

## The Math

In weighted least squares for a single protein, you minimize:

```
Σ_i  w_i * (y_i - x_i'β)²
```

If `w_i = w` (the same constant for every sample `i`), this becomes:

```
w * Σ_i (y_i - x_i'β)²
```

The constant `w` factors out of the optimization entirely. It doesn't change the location of the minimum, so:

- **β̂ is identical** to unweighted OLS
- **Residuals are identical**
- **σ² estimate is identical**
- **The QR decomposition is identical** (because `lm.wfit` absorbs weights into `sqrt(w) * X` and `sqrt(w) * y` — if `w` is constant, this is just a uniform rescaling that cancels)

So a per-protein constant weight literally does nothing in `lm.series`. It doesn't change the coefficients, doesn't change the standard errors, doesn't change the residual variance. It's a no-op.

And downstream in `eBayes`, the moderation works on `fit$sigma` and `fit$df.residual` — both unchanged — so the empirical Bayes step is also unaffected.

---

## Where Weights Actually Help

You need the weight to **vary across samples within a protein** for it to have any effect. This happens naturally in proteomics when you think about it correctly:

**The number of peptides quantified is usually not constant across samples.** Due to missing values at the peptide level, protein A might have 5 peptides quantified in sample 1, but only 2 in sample 3 (because 3 peptides had missing values in that sample). If you summarized by averaging, the protein-level mean in sample 1 is based on 5 peptides and is more precise than the mean in sample 3 based on 2 peptides.

So the correct weight is a **matrix**, not a vector:

```r
# weights[protein, sample] = number of non-missing peptides
#                             contributing to that protein's
#                             summary in that sample
weights <- matrix(..., nrow = n_proteins, ncol = n_samples)
fit <- lmFit(protein_summaries, design, weights = weights)
```

This gives sample 1's observation more influence than sample 3's for protein A, which is exactly what you want.

---

## What the Weight Should Actually Represent

If you're averaging `n` peptide log-intensities, each with variance `σ²_pep`, the variance of the mean is `σ²_pep / n`. The WLS weight should be proportional to precision (inverse variance), so:

```r
weight[protein, sample] = n_peptides[protein, sample]
```

This is correct under the assumption that peptides contribute independent, equally-variable measurements — which is a simplification (peptides vary in ionization efficiency, some are more noisy), but a reasonable first-order approximation.

More sophisticated weights could incorporate:

- **Per-peptide variance estimates** — if some peptides are consistently noisy, downweight them during summarization and propagate the resulting precision into the protein-level weight
- **maxLFQ-style weights** — maxLFQ doesn't simply average, so the effective "number of peptides" contributing to each protein-sample estimate is less straightforward, but the software can output uncertainty estimates

---

## The Case Where Peptide Counts Truly Are Constant

If you have no missingness at the peptide level (every protein has the same peptides quantified in every sample), then the weight degenerates to a constant per protein and is useless for `lmFit`. In that scenario, the only place the number of peptides *could* help is if you built it into a **prior** in `eBayes` — saying "I trust the variance estimate for proteins with many peptides more than for those with few." But standard `eBayes` doesn't have a mechanism for per-gene prior weights on the variance (it uses `df.residual` as the precision of each gene's variance estimate, which depends only on sample size, not peptide count).

---

## Weight Types in limma

| Weight type | What it does in limma | Useful? |
|---|---|---|
| Constant per protein (same for all samples) | Factors out of WLS, changes nothing | No |
| Varies per protein × sample (e.g., peptide count per sample) | Downweights imprecise observations | Yes |
| Array weights (same for all proteins in a sample) | Downweights noisy arrays | Yes, different use case |

The key insight: **weights only do something when they create contrast between observations within the same model** (i.e., across samples for a given protein). A uniform scaling per row is invisible to OLS.
