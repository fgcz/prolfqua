# How voom Integrates with limma

## TL;DR

Voom is a **preprocessing step** that estimates per-observation precision weights from a mean-variance trend, which are then passed as weights to limma's `lmFit` for weighted least squares. It is conceptually similar to a single iteration of IRLS in a GLM, but operates entirely within the linear model framework.

---

## What voom Does (Step by Step)

### Step 1: Transform Counts to log-CPM

Raw counts are converted to $\log_2$ counts-per-million:

```r
y <- t(log2(t(counts + 0.5) / (lib.size + 1) * 1e6))
```

The `+0.5` avoids $\log(0)$. This is the response variable that will ultimately be modelled by limma — voom transforms the data so we're working in continuous log-space, not in count space.

### Step 2: Fit a Preliminary Linear Model

Voom calls `lmFit(y, design)` on the log-CPM values. This gives preliminary estimates of the coefficients and, crucially, the **residual standard deviations** `fit$sigma` for each gene. This is an ordinary (or weighted) least squares fit — nothing fancy yet.

### Step 3: Estimate the Mean-Variance Trend

This is the core insight. Voom computes:

- $s_x$ = the average $\log_2$ count size per gene (gene-wise mean abundance)
- $s_y$ = $\sqrt{\hat{\sigma}_g}$ — the **square root of the residual standard deviation**

Note: `fit$sigma` in limma is the residual **standard deviation** $\sigma_g$, not the variance. So $s_y = \sigma_g^{1/2}$, which is equivalently the **fourth root of variance**: $s_y = (\sigma_g^2)^{1/4}$.

It then fits a **lowess smoother** through the $(s_x, s_y)$ scatter and wraps the result in a piecewise-linear interpolating function:

```r
sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e6)
sy <- sqrt(fit$sigma)
l  <- lowess(sx, sy, f = span)
f  <- approxfun(l, rule = 2, ties = list("ordered", mean))
```

The function $f$ maps any log-abundance value to the predicted $\sqrt{\text{sd}}$ at that abundance level. The `rule=2` means it extrapolates using the nearest value for points outside the range of the trend. This captures the fact that in RNA-seq, variance is not constant across genes — low-count genes are noisier than high-count genes, even on the log scale.

**Why smooth on the $\sqrt{\text{sd}}$ scale?** The distribution of residual standard deviations across genes is right-skewed. Taking the square root compresses that skew, making the lowess fit better behaved and more stable.

### Step 4: Compute Per-Observation Weights

Here's where it gets **observation-level**. Voom computes the fitted values from the preliminary model for every gene $\times$ sample combination:

```r
fitted.values   <- fit$coefficients %*% t(fit$design)
fitted.cpm      <- 2^fitted.values
fitted.count    <- 1e-6 * t(t(fitted.cpm) * (lib.size + 1))
fitted.logcount <- log2(fitted.count)
```

Then for each individual observation (gene $g$, sample $j$), it evaluates the lowess interpolation function $f$ at that observation's fitted log-count:

```r
w <- 1 / f(fitted.logcount)^4
```

#### Why `^4`?

The exponent traces back to the scale on which the lowess was fit:

$$f(\cdot) \;\text{returns}\; \sigma^{1/2} \;=\; (\sigma^2)^{1/4}$$

To recover the **variance** $\sigma^2$ (which is what we need for precision weights), we raise to the 4th power, then invert:

| Expression | What it is |
|-----------|-----------|
| $f(\hat{y}_{gj})$ | $\sigma_{gj}^{1/2}$ (fourth root of variance) |
| $f(\hat{y}_{gj})^2$ | $\sigma_{gj}$ (standard deviation) |
| $f(\hat{y}_{gj})^4$ | $\sigma_{gj}^2$ (variance) |
| $1 / f(\hat{y}_{gj})^4$ | $1/\sigma_{gj}^2$ (precision weight) |

So the final weight is:

$$w_{gj} = \frac{1}{\hat{\sigma}^2_{gj}}$$

> **Key point:** The weight for gene $g$ in sample $j$ depends on the *fitted count* for that specific observation, not just on the gene's average expression. A gene that is lowly expressed in one condition but highly expressed in another will receive different weights for samples in each condition.

### Step 5: Output

Voom returns an `EList` containing:

| Component     | Description                          |
|---------------|--------------------------------------|
| `$E`          | log-CPM expression matrix            |
| `$weights`    | Precision weight matrix (genes × samples) |
| `$design`     | Design matrix                        |

This `EList` object is what you pass directly to `lmFit`.

---

## How the Weights Enter the Linear Model

When you call `lmFit(v, design)`, limma detects the observation-level weights matrix and takes the **slow path** through `lm.series`. For each gene, it calls R's `lm.wfit(design, y_g, w_g)`, which solves the weighted least squares problem:

$$\hat{\beta}_g = (X^T W_g X)^{-1} \, X^T W_g \, y_g$$

where $W_g = \text{diag}(w_{g1}, \ldots, w_{gn})$ is the diagonal weight matrix for gene $g$.

Standard WLS: observations with higher predicted variance get downweighted.

### Why This Forces the Slow Path

Because the weights vary across **both genes and samples** (it's a full matrix, not a vector), the QR decomposition of $X^T W_g X$ is **different for every gene**. This prevents limma from using the fast single-QR path (which assumes a shared $(X^T X)^{-1}$ across all genes). Each gene effectively has its own design matrix $W_g^{1/2} X$.

### Covariance of the Coefficients

The covariance of the estimated coefficients for gene $g$ under WLS is:

$$\text{Cov}(\hat{\beta}_g) = \sigma^2_g \, (X^T W_g X)^{-1}$$

where $\sigma^2_g$ is estimated from the weighted residuals. After `eBayes()`, the moderated posterior variance $\tilde{\sigma}^2_g$ replaces $\sigma^2_g$.

---

## Relationship to GLMs

The connection is real, but voom is deliberately **not** a GLM. Here's the comparison:

### What a GLM Does (e.g., edgeR / DESeq2)

A Poisson or negative binomial GLM models counts directly with a log link. The mean-variance relationship is baked into the distributional assumption:

- Poisson: $\text{Var}(Y) = \mu$
- Negative binomial: $\text{Var}(Y) = \mu + \phi \mu^2$

Inference uses **iteratively reweighted least squares (IRLS)**, where at each iteration the working weights are derived from the current fitted values and the variance function.

### What Voom Does

Voom performs something conceptually similar but **stops after one iteration**:

1. **Transform to log scale** — analogous to the log-link step
2. **Fit a preliminary model** — one round of OLS
3. **Estimate variance as a smooth function of the fitted mean** — analogous to the "variance function," but estimated nonparametrically via lowess across all genes rather than assumed from a parametric family
4. **Compute weights** $= 1 / \text{variance}$ — the "working weights" step
5. **Fit WLS with those fixed weights** — one round of WLS, **not iterated**

### Key Differences

| Aspect | GLM (edgeR/DESeq2) | Voom + limma |
|--------|-------------------|--------------|
| **Response** | Raw counts | log-CPM (continuous) |
| **Variance function** | Parametric (NB family) | Nonparametric (lowess trend) |
| **Estimation** | IRLS (iterative) | Single iteration |
| **Per-gene dispersion** | Estimated + shrunk per gene | Implicitly captured via trend |
| **Downstream tools** | GLM-specific | Full limma ecosystem |

### Why One Iteration Is Sufficient

Law et al. (2014) argued that one iteration suffices because the $\log$-CPM transformation already stabilises the mean-variance relationship enough that the initial OLS estimates are reasonable starting points. The lowess smoothing across thousands of genes also provides a robust estimate of the trend.

---

## Architectural Advantages

The beauty of voom is that it converts the "RNA-seq count modelling problem" into a "weighted linear model problem." This means you inherit the entire limma ecosystem for free:

- Flexible design matrices (any model you can write as $X\beta$)
- Blocking and duplicate correlation (`gls.series`)
- Arbitrary contrasts
- Empirical Bayes moderation (`eBayes`)
- `treat` for fold-change thresholds
- Gene set testing (`camera`, `roast`, etc.)

GLM-based approaches (edgeR, DESeq2) need to re-implement all of that machinery within the GLM framework.

### Tradeoffs

Voom's variance estimates are approximate (one iteration, smoothed across genes) whereas edgeR/DESeq2 estimate per-gene dispersions with shrinkage. In practice, the two approaches perform very similarly for well-expressed genes with reasonable library sizes. Voom can struggle with **very sparse data** (many zeros), which is why `voomLmFit` in edgeR was developed as a more robust variant.

---

## Implementation Summary

If implementing voom + limma, the integration point is clean:

1. **Voom produces** a weight matrix $W$ of shape (genes $\times$ samples)
2. **`lmFit` solves** per-gene WLS using those weights:

$$\hat{\beta}_g = (X^T W_g X)^{-1} X^T W_g y_g$$

3. **Residual variance** is estimated from weighted residuals:

$$\hat{\sigma}^2_g = \frac{(y_g - X\hat{\beta}_g)^T W_g (y_g - X\hat{\beta}_g)}{n - p}$$

4. **Everything downstream** (`eBayes`, contrasts, etc.) is unchanged — it uses the WLS residual variances and the per-gene $(X^T W_g X)^{-1}$ as the unscaled covariance.

---

## Extension to Proteomics (Already log2-Transformed Data)

### The Procedure

For proteomics data that is already log2-transformed, the procedure simplifies because there is no count-to-logCPM transformation and no library size offset. Limma provides `vooma()` for this purpose, but the logic is straightforward to understand:

1. **`lmFit(y, design)`** — ordinary least squares on log2 intensities, **no weights**
2. **Lowess on `(fit$Amean, sqrt(fit$sigma))`** — estimate the mean-variance trend
3. **Evaluate the trend at each observation's fitted value** → derive precision weights
4. **`lmFit(y, design, weights=w)`** — WLS refit with the estimated weights

### How Fitted Values Are Computed (No Library Size)

In RNA-seq voom, there is a back-transformation step that re-introduces a **per-sample library size**. `lib.size` is a vector of length $n$ (one value per sample), representing the total number of reads sequenced for that sample:

```r
lib.size <- colSums(counts)
# e.g., lib.size = c(25e6, 32e6, 18e6, 28e6, ...)
```

This is what makes RNA-seq voom weights unique per observation even within the same group:

```r
fitted.values  <- fit$coefficients %*% t(fit$design)     # logCPM — same within group
fitted.cpm     <- 2^fitted.values                          # CPM — same within group  
fitted.count   <- 1e-6 * t(t(fitted.cpm) * (lib.size + 1)) # COUNT — differs per sample!
fitted.logcount <- log2(fitted.count)                       # log-count — differs per sample
```

Even though two samples are in the same group (same fitted CPM), they have different total read depths, so their fitted *counts* differ, and therefore their weights differ. A sample sequenced to 50M reads gets a different weight than one at 20M reads, even for the same gene in the same group.

**For proteomics, none of this applies.** The data is already on a log2 intensity scale with no sample-specific offset. The fitted values come directly from the linear model:

```r
fitted.values <- fit$coefficients %*% t(fit$design)
```

That's it. For a two-group intercept + treatment design:

$$\hat{y}_{gj} = \begin{cases} \hat{\beta}_{g0} & \text{if sample } j \in \text{group 1} \\[4pt] \hat{\beta}_{g0} + \hat{\beta}_{g1} & \text{if sample } j \in \text{group 2} \end{cases}$$

You evaluate the lowess trend $f$ at these fitted values to get the predicted $\sqrt{\text{sd}}$, then:

$$w_{gj} = \frac{1}{f(\hat{y}_{gj})^4}$$

Since $\hat{y}_{gj}$ is identical for all samples in the same group, $w_{gj}$ is identical too.

### Complete Proteomics Procedure (Manual)

```r
# Step 1: OLS fit (no weights)
# vooma uses lm.fit directly (not lmFit), operating on transposed data
fit <- lm.fit(design, t(y))
mu  <- fit$fitted.values              # fitted values: samples × proteins

# Step 2: Compute residual variance from QR effects
# s2 is the residual variance per protein, computed from the
# orthogonalized residual effects (more numerically stable than RSS)
s2 <- colMeans(fit$effects[-(1:fit$rank), , drop = FALSE]^2)

# Step 3: Lowess trend on fourth-root-of-variance scale
sx <- rowMeans(y)                     # mean log2-intensity per protein
sy <- sqrt(sqrt(s2))                  # (variance)^(1/4) = (sd)^(1/2)
l  <- lowess(sx, sy, f = 0.5)
f  <- approxfun(l, rule = 2)         # piecewise-linear interpolation of the trend

# Step 4: Observation-level precision weights
# f returns (variance)^(1/4), so f^4 = variance, and 1/f^4 = precision
w <- 1 / f(t(mu))^4                  # t(mu) puts it back to proteins × samples
dim(w) <- dim(y)

# Step 5: WLS refit with voom-style weights
fit2 <- lmFit(y, design, weights = w)
fit2 <- eBayes(fit2)
```

### Why `vooma()` Uses `lm.fit` + `s2` Instead of `lmFit` + `fit$sigma`

If you read the `vooma()` source, you'll notice it computes residual variances via:

```r
fit <- lm.fit(design, t(y$E))
s2  <- colMeans(fit$effects[-(1:fit$rank), , drop = FALSE]^2)
sy  <- sqrt(sqrt(s2))
```

whereas the manual approach uses:

```r
fit <- lmFit(y, design)
sy  <- sqrt(fit$sigma)
```

These are **computationally identical**, not just mathematically equivalent. Here's why.

When `lmFit` calls `lm.series`, which calls `lm.fit`, the residual standard deviation `fit$sigma` is computed from exactly the same QR effects:

```r
# Inside lm.series (called by lmFit):
effects <- fit$effects
sigma   <- sqrt(colMeans(effects[-(1:fit$rank), , drop = FALSE]^2))
```

So `fit$sigma` from `lmFit` is $\sqrt{s^2}$, and the chain of equivalences is:

| Expression | Equals | What it is |
|-----------|--------|-----------|
| `s2` | $s^2$ | Residual variance per protein |
| `fit$sigma` | $\sqrt{s^2}$ | Residual standard deviation (same `s2`, just square-rooted) |
| `sqrt(sqrt(s2))` | $(s^2)^{1/4}$ | Fourth root of variance |
| `sqrt(fit$sigma)` | $(\sqrt{s^2})^{1/2} = (s^2)^{1/4}$ | Same fourth root of variance |

The difference between the two approaches is purely about **efficiency**, not numerical quality:

- `vooma()` uses `lm.fit` directly because it's computing a throwaway preliminary fit. It only needs fitted values and residual variances, then discards everything. Calling `lmFit` would construct the full `MArrayLM` object — gene annotation extraction via `getEAWP`, coefficient names, covariance structures, etc. — all overhead that gets thrown away immediately.
- The `lmFit` + `fit$sigma` approach is slightly less efficient but perfectly correct, and it has a practical advantage: `lmFit` can handle observation-level weight matrices (needed when incorporating external weights like precursor counts), whereas `lm.fit` cannot.

**Bottom line:** For custom implementations — especially when you need external weights — using `lmFit` + `sqrt(fit$sigma)` is the right choice. There is no numerical downside.

### Using `vooma()` Directly

Limma provides `vooma()` which wraps the above procedure in a single call:

```r
v <- vooma(y, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
```

There is also `voomaByGroup()`, which fits **separate mean-variance trends per group**. This is useful if different experimental conditions have different variance structures (e.g., treatment samples are systematically noisier than controls):

```r
v <- voomaByGroup(y, group = condition, design = design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
```

> **Important limitation:** Unlike `voom()`, `vooma()` does **not** accept a `weights` argument.
> If you need to combine voom-style weights with external weights (e.g., precursor counts),
> you must use the manual procedure — see the section on integrating external weights below.

### Weight Structure in Proteomics

For a simple two-group design with 4 samples per group, the weight vector for protein $g$ looks like:

$$w_g = (\underbrace{x, x, x, x}_{\text{group 1}}, \; \underbrace{y, y, y, y}_{\text{group 2}})$$

where $x$ and $y$ are determined by the fitted log2-intensity for protein $g$ in each group. The weights **do** vary across proteins — each protein sits at a different position on the mean-variance trend curve — but within a group, all samples for the same protein share the same weight.

With more complex designs (continuous covariates, interaction terms, batch effects as continuous variables), each sample could receive a unique fitted value and therefore a unique weight, even for proteomics.

### Comparison: RNA-seq voom vs. Proteomics voom

| Aspect | RNA-seq voom | Proteomics (vooma) |
|--------|-------------|-------------------|
| **Input** | Raw counts | log2 intensities |
| **log-CPM transform** | Yes (`log2(count + 0.5)`) | Not needed |
| **Library size offset** | Yes — per-sample $L_j$ | No |
| **Fitted values for trend** | $\log_2(\hat{\mu}_{gj} \cdot L_j)$ — unique per observation | $X\hat{\beta}_g$ — same within group (simple designs) |
| **Weight granularity** | Every cell in the matrix can be unique | At most $k$ unique values per protein ($k$ = number of distinct fitted values) |

---

## Integrating External Weights (e.g., Precursor Counts)

In proteomics, you may have additional precision information beyond the mean-variance relationship — for example, the number of precursors (peptides) contributing to each protein's intensity estimate. Proteins quantified from 20 precursors are more precisely estimated than those from 2. This information can be combined with voom-style weights.

### Why Multiplication Is Correct

The two weight sources capture **independent** precision information:

- **Voom weights** capture: "how noisy is a protein at this abundance level?" — the mean-variance relationship
- **Precursor weights** capture: "how well-estimated is this protein's intensity?" — more precursors = more peptide evidence = more precise

If the variance of observation $(g, j)$ is $\sigma^2_{gj}$ from the mean-variance trend, and the measurement precision is proportional to $n_{gj}$ (number of precursors), then the effective variance is:

$$\text{Var}_{\text{eff}}(y_{gj}) \propto \frac{\sigma^2_{gj}}{n_{gj}}$$

So the precision weight is:

$$w_{gj} = \frac{n_{gj}}{\sigma^2_{gj}} = n_{gj} \times w^{\text{voom}}_{gj}$$

Which is the product of the two weight sources.

### Where External Weights Enter

External weights should enter at **two** points:

1. **Into the preliminary fit** used to estimate the mean-variance trend. This way the trend estimation itself accounts for the fact that some observations are more precise than others.
2. **Into the final WLS fit**, by multiplying voom weights $\times$ external weights.

### Complete Procedure with External Weights

Since `vooma()` does not accept external weights, this must be done manually:

```r
# External weights: e.g., number of precursors per protein × sample
# Shape: same as expression matrix (proteins × samples)
# If per-protein only (not varying across samples), expand to a matrix:
#   precursor_weights <- matrix(n_precursors, nrow=nrow(y), ncol=ncol(y))

# Step 1: Initial fit incorporating external weights
# This ensures the mean-variance trend is estimated from weighted residuals
fit <- lmFit(y, design, weights = precursor_weights)

# Step 2: Lowess trend on fourth-root-of-variance scale
# fit$sigma from lmFit is the residual SD; sqrt gives (variance)^(1/4)
sx <- fit$Amean                      # = rowMeans(y), mean log2-intensity per protein
sy <- sqrt(fit$sigma)                # equivalent to sqrt(sqrt(s2)) = (variance)^(1/4)
l  <- lowess(sx, sy, f = 0.5)
f  <- approxfun(l, rule = 2)

# Step 3: Voom-style weights from the trend
fitted.values <- fit$coefficients %*% t(fit$design)
w_voom <- 1 / f(fitted.values)^4    # f^4 = variance, 1/f^4 = precision
dim(w_voom) <- dim(fitted.values)

# Step 4: Combined weights = product of both sources
w_combined <- w_voom * precursor_weights

# Step 5: Final WLS fit with combined weights
fit2 <- lmFit(y, design, weights = w_combined)
fit2 <- eBayes(fit2)
```

> **Note:** Here we use `lmFit` (not `lm.fit`) because we need it to handle the external
> weight matrix — `lm.fit` cannot accept observation-level weight matrices in the way limma
> structures them. The `fit$sigma` path is computationally identical to vooma's `sqrt(sqrt(s2))`
> — see "Why `vooma()` Uses `lm.fit` + `s2` Instead of `lmFit` + `fit$sigma`" above for the
> full derivation.

### Practical Considerations for Precursor Weights

Raw precursor counts can produce very extreme weight ratios. A protein with 1 precursor vs. one with 50 precursors would get a 50× difference in weight before the voom component is even applied. This can cause a small number of high-precursor proteins to dominate the fit. Common strategies to mitigate this:

| Strategy | Code | Effect |
|----------|------|--------|
| **Cap at a maximum** | `pmin(n_precursors, 20)` | Limits influence of very high-count proteins |
| **Square root** | `sqrt(n_precursors)` | Compresses the range (1→1, 4→2, 25→5, 50→7) |
| **Log transform** | `log(n_precursors + 1)` | Aggressive compression (1→0.7, 10→2.4, 50→3.9) |
| **Rank-based** | `rank(n_precursors) / max(rank(...))` | Uniform spread, removes outlier influence |

The right choice depends on whether you believe precision truly scales linearly with precursor count (use raw or capped) or whether the relationship saturates (use sqrt or log). In practice, `sqrt` is often a reasonable default — it acknowledges that going from 1 to 4 precursors matters more than going from 46 to 49.

### Note on Weight Dimensions

External weights can be:

- **A matrix** (proteins × samples): when precursor counts vary per protein *and* per sample (e.g., different peptides detected in different runs). This is the most common case in DIA or MaxQuant-style output.
- **A vector** (one per protein): when you have a single precursor count per protein, expanded to a matrix via `matrix(n_precursors, nrow=nrow(y), ncol=ncol(y))`. This assumes the same measurement quality across all samples for a given protein.

When the external weights vary per sample, the combined weight matrix gains additional granularity beyond what voom alone provides — now every cell can be unique even for simple group designs.

---

## References

- Law, C.W., Chen, Y., Shi, W. & Smyth, G.K. (2014). Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology* 15, R29.
- Smyth, G.K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Statistical Applications in Genetics and Molecular Biology* 3, Article 3.
