# Decorate LFQData with Methods for transforming Intensities

Transformations of the response of an
[LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object:
log2, robust scaling (`robscale()`), centring to a reference subset,
subtraction of sample offsets and any function applied to the intensity
array or matrix. Obtained with `LFQData$get_Transformer()`; methods
return the transformer, the result is in `$lfq`.

## Value

An R6 class generator.

## Centring and contrasts

Centring is an explicit transformation step, never applied by a model.
`robscale()` and `robscale_subset()` centre each sample at its median
(and scale by its MAD); `preserve_mean = TRUE` adds the overall mean
back, so abundances stay on the original level, `FALSE` centres them at
zero. `center_to_reference()` and `subtract_sample_offsets()` shift each
sample by one offset. A per-sample shift moves every feature of that
sample by the same amount: a contrast between groups changes by the
difference of the mean offsets of the groups, and a shift common to all
samples (such as `preserve_mean`) changes the average abundance
(`avgAbd`) but not the contrast estimates, standard errors or p-values.

## Public fields

- `lfq`:

  LFQData object

## Methods

### Public methods

- [`LFQDataTransformer$new()`](#method-LFQDataTransformer-initialize)

- [`LFQDataTransformer$log2()`](#method-LFQDataTransformer-log2)

- [`LFQDataTransformer$get_scales()`](#method-LFQDataTransformer-get_scales)

- [`LFQDataTransformer$robscale()`](#method-LFQDataTransformer-robscale)

- [`LFQDataTransformer$robscale_subset()`](#method-LFQDataTransformer-robscale_subset)

- [`LFQDataTransformer$center_to_reference()`](#method-LFQDataTransformer-center_to_reference)

- [`LFQDataTransformer$subtract_sample_offsets()`](#method-LFQDataTransformer-subtract_sample_offsets)

- [`LFQDataTransformer$intensity_array()`](#method-LFQDataTransformer-intensity_array)

- [`LFQDataTransformer$intensity_matrix()`](#method-LFQDataTransformer-intensity_matrix)

- [`LFQDataTransformer$clone()`](#method-LFQDataTransformer-clone)

------------------------------------------------------------------------

### `LFQDataTransformer$new()`

initialize

#### Usage

    LFQDataTransformer$new(lfqdata)

#### Arguments

- `lfqdata`:

  LFQData object to transform

------------------------------------------------------------------------

### `LFQDataTransformer$log2()`

log2 transform data

#### Usage

    LFQDataTransformer$log2(force = FALSE)

#### Arguments

- `force`:

  if `FALSE`, data already log2 transformed will not be transformed a
  second time. `TRUE` forces log transformation

#### Returns

LFQDataTransformer

------------------------------------------------------------------------

### `LFQDataTransformer$get_scales()`

get mean and variance and standard deviation in each sample

#### Usage

    LFQDataTransformer$get_scales()

#### Returns

list with means and mads

------------------------------------------------------------------------

### `LFQDataTransformer$robscale()`

robust scale data

#### Usage

    LFQDataTransformer$robscale(
      preserve_mean = TRUE,
      colname = "transformedIntensity"
    )

#### Arguments

- `preserve_mean`:

  should original mean value be preserved TRUE, if FALSE then center at
  zero; changes the average abundance, not the contrasts (see section
  "Centring and contrasts")

- `colname`:

  new name of transformed column

#### Returns

LFQDataTransformer (self)

------------------------------------------------------------------------

### `LFQDataTransformer$robscale_subset()`

log2 transform and robust scale data based on subset

#### Usage

    LFQDataTransformer$robscale_subset(
      lfqsubset,
      preserve_mean = TRUE,
      colname = "transformed_abundance"
    )

#### Arguments

- `lfqsubset`:

  LFQData subset to use for normalization

- `preserve_mean`:

  should original mean value be preserved TRUE, if FALSE then center at
  zero

- `colname`:

  how to name the transformed intensities, default transformedIntensity

#### Returns

LFQDataTransformer (self)

------------------------------------------------------------------------

### `LFQDataTransformer$center_to_reference()`

center data to a reference subset

#### Usage

    LFQDataTransformer$center_to_reference(lfqsubset)

#### Arguments

- `lfqsubset`:

  LFQData subset to use as reference

#### Returns

LFQDataTransformer (self)

------------------------------------------------------------------------

### `LFQDataTransformer$subtract_sample_offsets()`

Subtract caller-supplied offsets from the active log2 response. The
caller supplies log2 data and chooses the global level through the
offsets. Missing intensities remain missing. The result becomes the
active response in a new column named `<response>_offset_adjusted`, with
a numeric suffix if the name already exists.

#### Usage

    LFQDataTransformer$subtract_sample_offsets(offsets)

#### Arguments

- `offsets`:

  Named numeric vector with one finite log2 offset per sample. Names
  must match the configured sample IDs exactly and may be in any order.

#### Returns

LFQDataTransformer (self, invisibly)

------------------------------------------------------------------------

### `LFQDataTransformer$intensity_array()`

Transforms intensities

#### Usage

    LFQDataTransformer$intensity_array(
      .func = log2,
      force = FALSE,
      transformed = TRUE
    )

#### Arguments

- `.func`:

  transformation function working with arrays e.g. log2, log10, asinh
  etc.

- `force`:

  transformation on already transformed data.

- `transformed`:

  whether the result counts as transformed; set `FALSE` for a
  back-transformation such as `exp`

#### Returns

LFQDataTransformer (self)

------------------------------------------------------------------------

### `LFQDataTransformer$intensity_matrix()`

pass a function which works with matrices, e.g., vsn::justvsn

#### Usage

    LFQDataTransformer$intensity_matrix(.func = .robust_scale, force = FALSE)

#### Arguments

- `.func`:

  any function taking a matrix and returning a matrix (columns sample,
  rows feature, e.g.
  [`base::scale()`](https://rdrr.io/r/base/scale.html)); default: robust
  scaling, which centres each sample on its median and scales it by its
  MAD, as `$robscale()`

- `force`:

  transformation on data already transformed

#### Returns

LFQDataTransformer (self)

------------------------------------------------------------------------

### `LFQDataTransformer$clone()`

The objects of this class are cloneable with this method.

#### Usage

    LFQDataTransformer$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config(n_proteins = 50)
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 225 rows with NA response for feature and sample combinations absent from the input.
lfqdata <- LFQData$new(istar$data, istar$config)

lfqcopy <- lfqdata$get_copy()
lfqTrans <- lfqcopy$get_Transformer()

x <- lfqTrans$intensity_array(log2)
#> Transformed response column added: log2_abundance

x$lfq$is_transformed()
#> [1] TRUE
x <- x$intensity_matrix()
#> Warning: The data are already transformed; set `force = TRUE` to log2-transform them again.
plotter <- x$lfq$get_Plotter()
plotter$intensity_distribution_density()


# transform by asinh root and scale

lfqcopy <- lfqdata$get_copy()
lfqTrans <- lfqcopy$get_Transformer()
x <- lfqTrans$intensity_array(asinh)
#> Transformed response column added: asinh_abundance
mads1 <- mean(x$get_scales()$mads)
x <- lfqTrans$intensity_matrix(force = TRUE)
mads2 <- mean(x$get_scales()$mads)
stopifnot(abs(mads1 - mads2) < 1e-8)


stopifnot(abs(mean(x$get_scales()$medians)) < 1e-8)
lfqcopy <- lfqdata$get_copy()
lfqTrans <- lfqcopy$get_Transformer()
lfqTrans$log2()
#> Transformed response column added: log2_abundance
before <- lfqTrans$get_scales()
lfqTrans$robscale()
after <- lfqTrans$get_scales()
stopifnot(abs(mean(before$medians) - mean(after$medians)) < 1e-8)
stopifnot(abs(mean(before$mads) - mean(after$mads)) < 1e-8)

# normalize data using vsn
lfqcopy <- lfqdata$get_copy()
lfqTrans <- lfqcopy$get_Transformer()
lfqTransCheck <- lfqcopy$get_Transformer()

lfqTransCheck$log2()
#> Transformed response column added: log2_abundance
lfqTransCheck$get_scales()
#> $medians
#>     A_V1     A_V2     A_V3     A_V4     B_V1     B_V2     B_V3     B_V4 
#> 4.463740 4.450835 4.444549 4.463231 4.537348 4.568077 4.495236 4.519219 
#>  Ctrl_V1  Ctrl_V2  Ctrl_V3  Ctrl_V4 
#> 4.500299 4.453443 4.491344 4.464994 
#> 
#> $mads
#>      A_V1      A_V2      A_V3      A_V4      B_V1      B_V2      B_V3      B_V4 
#> 0.3515227 0.3699336 0.3679643 0.3853720 0.3884622 0.3747087 0.3602410 0.3672785 
#>   Ctrl_V1   Ctrl_V2   Ctrl_V3   Ctrl_V4 
#> 0.3507487 0.3606693 0.3632381 0.3691903 
#> 
lfqTransCheck$lfq$get_Plotter()$intensity_distribution_density()


if (requireNamespace("vsn", quietly = TRUE)) {
 res <- lfqTrans$intensity_matrix( .func = vsn::justvsn)
 res$lfq$get_Plotter()$intensity_distribution_density()
 res$get_scales()
}
#> $medians
#>     A_V1     A_V2     A_V3     A_V4     B_V1     B_V2     B_V3     B_V4 
#> 5.010572 5.002975 4.994914 5.006013 5.051370 5.074292 5.033370 5.044391 
#>  Ctrl_V1  Ctrl_V2  Ctrl_V3  Ctrl_V4 
#> 5.043057 5.010404 5.033362 5.015886 
#> 
#> $mads
#>      A_V1      A_V2      A_V3      A_V4      B_V1      B_V2      B_V3      B_V4 
#> 0.2443545 0.2638516 0.2494166 0.2556915 0.2710033 0.2638501 0.2316761 0.2437389 
#>   Ctrl_V1   Ctrl_V2   Ctrl_V3   Ctrl_V4 
#> 0.2738565 0.2415127 0.2514121 0.2512947 
#> 
if (requireNamespace("preprocessCore", quietly = TRUE)) {
quant <- function(y){
 ynorm <- preprocessCore::normalize.quantiles(y)
 rownames(ynorm) <- rownames(y)
 colnames(ynorm) <- colnames(y)
 return(ynorm)
}
 res <- lfqTrans$intensity_matrix( .func = quant)
 res$lfq$get_Plotter()$intensity_distribution_density()
}
#> Warning: The data are already transformed; set `force = TRUE` to log2-transform them again.



#subset scaling

istar2 <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.
lfqdata2 <- LFQData$new(istar2$data, istar2$config)
lfqdata2 <- lfqdata2$get_Transformer()$intensity_array(log2)$lfq
#> Transformed response column added: log2_abundance
head(lfqdata2$hierarchy())
#> # A tibble: 6 × 1
#>   protein_Id 
#>   <chr>      
#> 1 0EfVhX~0087
#> 2 7cbcrd~5725
#> 3 9VUkAq~4703
#> 4 BEJI92~5282
#> 5 CGzoYe~2147
#> 6 DoWup2~5896
internal <- lfqdata2$get_subset(head(lfqdata2$hierarchy()))
internal$hierarchy()
#> # A tibble: 6 × 1
#>   protein_Id 
#>   <chr>      
#> 1 0EfVhX~0087
#> 2 7cbcrd~5725
#> 3 9VUkAq~4703
#> 4 BEJI92~5282
#> 5 CGzoYe~2147
#> 6 DoWup2~5896
tr <- lfqdata2$get_Transformer()
tr$center_to_reference(internal)
pl <- tr$lfq$get_Plotter()
pl$intensity_distribution_density()

lfqdata2$get_Plotter()$intensity_distribution_density()

robscale <- lfqdata2$get_Transformer()$robscale_subset(internal)$lfq
robscale$get_Plotter()$intensity_distribution_density()


# Apply a one-unit log2 offset to the first sample
tr <- lfqdata$get_Transformer()$log2()
#> Transformed response column added: log2_abundance
samples <- unique(tr$lfq$data_long()[[tr$lfq$sample_name()]])
offsets <- setNames(numeric(length(samples)), samples)
offsets[1] <- 1
tr$subtract_sample_offsets(offsets)
```
