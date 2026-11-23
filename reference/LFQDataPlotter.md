# LFQDataPlotter —- Create various visualization of the LFQdata

LFQDataPlotter —- Create various visualization of the LFQdata

## Value

An R6 class generator.

## See also

Other LFQData:
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md),
[`AggregateMedpolish`](https://fgcz.github.io/prolfqua/reference/AggregateMedpolish.md),
[`AggregateRlm`](https://fgcz.github.io/prolfqua/reference/AggregateRlm.md),
[`AggregateTopN`](https://fgcz.github.io/prolfqua/reference/AggregateTopN.md),
[`AggregatorBase`](https://fgcz.github.io/prolfqua/reference/AggregatorBase.md),
[`LFQData`](https://fgcz.github.io/prolfqua/reference/LFQData.md),
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://fgcz.github.io/prolfqua/reference/LFQDataSummariser.md),
[`lfqdata_from_wide()`](https://fgcz.github.io/prolfqua/reference/lfqdata_from_wide.md),
[`lfqdata_to_summarized_experiment()`](https://fgcz.github.io/prolfqua/reference/lfqdata_to_summarized_experiment.md),
[`plot.LFQData()`](https://fgcz.github.io/prolfqua/reference/plot.LFQData.md)

## Public fields

- `lfq`:

  LFQData object

- `prefix`:

  prefix to figure names when writing, e.g. protein\_

## Methods

### Public methods

- [`LFQDataPlotter$new()`](#method-LFQDataPlotter-initialize)

- [`LFQDataPlotter$raster()`](#method-LFQDataPlotter-raster)

- [`LFQDataPlotter$heatmap()`](#method-LFQDataPlotter-heatmap)

- [`LFQDataPlotter$heatmap_cor()`](#method-LFQDataPlotter-heatmap_cor)

- [`LFQDataPlotter$pca()`](#method-LFQDataPlotter-pca)

- [`LFQDataPlotter$pca_plotly()`](#method-LFQDataPlotter-pca_plotly)

- [`LFQDataPlotter$boxplots()`](#method-LFQDataPlotter-boxplots)

- [`LFQDataPlotter$missingness_histogram()`](#method-LFQDataPlotter-missingness_histogram)

- [`LFQDataPlotter$missigness_histogram()`](#method-LFQDataPlotter-missigness_histogram)

- [`LFQDataPlotter$na_heatmap()`](#method-LFQDataPlotter-na_heatmap)

- [`LFQDataPlotter$intensity_distribution_density()`](#method-LFQDataPlotter-intensity_distribution_density)

- [`LFQDataPlotter$intensity_distribution_violin()`](#method-LFQDataPlotter-intensity_distribution_violin)

- [`LFQDataPlotter$pairs_smooth()`](#method-LFQDataPlotter-pairs_smooth)

- [`LFQDataPlotter$sample_correlation()`](#method-LFQDataPlotter-sample_correlation)

- [`LFQDataPlotter$upset_missing()`](#method-LFQDataPlotter-upset_missing)

- [`LFQDataPlotter$write_boxplots()`](#method-LFQDataPlotter-write_boxplots)

- [`LFQDataPlotter$write_pdf()`](#method-LFQDataPlotter-write_pdf)

- [`LFQDataPlotter$clone()`](#method-LFQDataPlotter-clone)

------------------------------------------------------------------------

### `LFQDataPlotter$new()`

create LFQDataPlotter

#### Usage

    LFQDataPlotter$new(lfqdata, prefix = "ms_")

#### Arguments

- `lfqdata`:

  [LFQData](https://fgcz.github.io/prolfqua/reference/LFQData.md) object
  with the data and its configuration

- `prefix`:

  will be prepended to outputs written

------------------------------------------------------------------------

### `LFQDataPlotter$raster()`

plot intensities in raster

#### Usage

    LFQDataPlotter$raster(
      arrange = c("mean", "var"),
      not_na = FALSE,
      rownames = FALSE,
      max_rownames_chars = 60,
      max_sample_label_chars = 20
    )

#### Arguments

- `arrange`:

  arrange by either mean or var

- `not_na`:

  TRUE arrange by number of NA's, FALSE by arrange by intensity

- `rownames`:

  show rownames (default FALSE - do not show.)

- `max_rownames_chars`:

  maximum displayed row label length

- `max_sample_label_chars`:

  maximum displayed sample label length. Labels keep their suffix
  because sample prefixes are often shared.

#### Returns

ComplexHeatmap::Heatmap

------------------------------------------------------------------------

### `LFQDataPlotter$heatmap()`

heatmap of intensities - columns are samples, rows are proteins or
peptides.

The abundances of each protein (row) are z-scored. Afterward, the mean
abundance for each protein is zero, and the standard variation is one.
z-scoring allows to compare (cluster) the proteins according to the
difference in the expression in the samples. Without the z-scoring, the
proteins would group according to their abundance, e.g., high abundant
proteins would be one cluster.

By default only the `n_most_variable` most variable features are shown
and the rows are clustered. Row clustering uses
[`stats::hclust`](https://rdrr.io/r/stats/hclust.html), which errors
above 65536 features, so the rows are ranked by their variability
statistic (CV for untransformed data, sd for transformed data; see
[`LFQDataStats`](https://fgcz.github.io/prolfqua/reference/LFQDataStats.md))
and the most variable are kept. This keeps the heatmap feasible and
legible for peptide-list / entrapment searches with tens of thousands of
features.

With `row_order`, the heatmap shows exactly these features in this
order, without clustering. Use it to compare two states of the same data
row by row, e.g. before and after imputation: imputation changes the
per-feature sd, so the most variable set and the clustering would differ
between the two heatmaps.

#### Usage

    LFQDataPlotter$heatmap(
      na_fraction = 1,
      rownames = FALSE,
      max_rownames_chars = 60,
      max_sample_label_chars = 20,
      n_most_variable = 1000,
      row_order = NULL,
      top_n
    )

#### Arguments

- `na_fraction`:

  maximum fraction of missing values per row. The default keeps every
  row that can be meaningfully z-scored.

- `rownames`:

  show rownames (default FALSE - do not show.)

- `max_rownames_chars`:

  maximum displayed row label length

- `max_sample_label_chars`:

  maximum displayed sample label length. Labels keep their suffix
  because sample prefixes are often shared.

- `n_most_variable`:

  keep the `n_most_variable` most variable features (default 1000);
  `NULL` or `Inf` keeps all features. Not together with `row_order`.

- `row_order`:

  feature IDs (rownames of `LFQData$data_wide()`) to show, in this
  order, without row clustering; default `NULL` selects and clusters the
  most variable features.

- `top_n`:

  deprecated, use `n_most_variable`

#### Returns

ComplexHeatmap::Heatmap

------------------------------------------------------------------------

### `LFQDataPlotter$heatmap_cor()`

heatmap of sample correlations.

The Spearman correlation among all samples is computed. Then the
euclidean distance is used to compute the distances.

#### Usage

    LFQDataPlotter$heatmap_cor(max_sample_label_chars = 20, method = "pearson")

#### Arguments

- `max_sample_label_chars`:

  maximum displayed sample label length. Labels keep their suffix
  because sample prefixes are often shared.

- `method`:

  correlation method: `"pearson"` (default), `"spearman"` or `"kendall"`
  for [`stats::cor()`](https://rdrr.io/r/stats/cor.html) on pairwise
  complete observations, or a function of the matrix returning a
  sample-by-sample correlation matrix

#### Returns

ComplexHeatmap::Heatmap

------------------------------------------------------------------------

### `LFQDataPlotter$pca()`

PCA plot

A PCA is applied and the first and second principal component are shown.
Features with missing values are removed. To keep all features, impute
first, e.g. with
[`AggregateLimpa`](https://fgcz.github.io/prolfqua/reference/AggregateLimpa.md)
in `impute_only` mode.

#### Usage

    LFQDataPlotter$pca(PC = c(1, 2), add_txt = TRUE, nudge = 0.1)

#### Arguments

- `PC`:

  default c(1,2) - first and second principal component

- `add_txt`:

  show sample names

- `nudge`:

  default 0.1 nudge point labels

#### Returns

ggplot

------------------------------------------------------------------------

### `LFQDataPlotter$pca_plotly()`

pca plot

#### Usage

    LFQDataPlotter$pca_plotly(PC = c(1, 2), add_txt = FALSE)

#### Arguments

- `PC`:

  default c(1,2) - first and second principal component

- `add_txt`:

  show sample names

#### Returns

plotly

------------------------------------------------------------------------

### `LFQDataPlotter$boxplots()`

boxplots for all proteins

#### Usage

    LFQDataPlotter$boxplots(facet = TRUE)

#### Arguments

- `facet`:

  enable facet wrap if hierarchy_depth less then hierarchy length.

#### Returns

tibble with column boxplots containing ggplot objects

------------------------------------------------------------------------

### `LFQDataPlotter$missingness_histogram()`

histogram of intensities given number of missing in conditions

#### Usage

    LFQDataPlotter$missingness_histogram()

#### Returns

ggplot

------------------------------------------------------------------------

### `LFQDataPlotter$missigness_histogram()`

Deprecated; use `$missingness_histogram()`.

#### Usage

    LFQDataPlotter$missigness_histogram()

#### Returns

ggplot

------------------------------------------------------------------------

### `LFQDataPlotter$na_heatmap()`

heatmap of features with missing values

#### Usage

    LFQDataPlotter$na_heatmap()

#### Returns

ComplexHeatmap::Heatmap

------------------------------------------------------------------------

### `LFQDataPlotter$intensity_distribution_density()`

density distribution of intensities

#### Usage

    LFQDataPlotter$intensity_distribution_density(
      legend = NA,
      max_legend_samples = 16
    )

#### Arguments

- `legend`:

  show legend TRUE, FALSE do not show, NA selects automatically based on
  sample count.

- `max_legend_samples`:

  maximum number of samples for automatic legend display.

#### Returns

ggplot

------------------------------------------------------------------------

### `LFQDataPlotter$intensity_distribution_violin()`

Violinplot showing distribution of intensities in all samples

#### Usage

    LFQDataPlotter$intensity_distribution_violin()

#### Returns

ggplot

------------------------------------------------------------------------

### `LFQDataPlotter$pairs_smooth()`

pairsplot of intensities

#### Usage

    LFQDataPlotter$pairs_smooth(max = 10)

#### Arguments

- `max`:

  maximal number of samples to show

#### Returns

NULL

------------------------------------------------------------------------

### `LFQDataPlotter$sample_correlation()`

plot of sample correlations

#### Usage

    LFQDataPlotter$sample_correlation(method = "pearson")

#### Arguments

- `method`:

  correlation method: `"pearson"` (default), `"spearman"` or `"kendall"`
  for [`stats::cor()`](https://rdrr.io/r/stats/cor.html) on pairwise
  complete observations, or a function of the matrix returning a
  sample-by-sample correlation matrix

#### Returns

NULL

------------------------------------------------------------------------

### `LFQDataPlotter$upset_missing()`

upset plot based on presence absence information

#### Usage

    LFQDataPlotter$upset_missing()

#### Returns

plot

------------------------------------------------------------------------

### `LFQDataPlotter$write_boxplots()`

write boxplots to file

#### Usage

    LFQDataPlotter$write_boxplots(path_qc, filename = NULL, width = 6, height = 6)

#### Arguments

- `path_qc`:

  path to write to

- `filename`:

  file name; one ending in `.pdf` is used as given, otherwise
  `boxplot.pdf` is appended (default: the plotter prefix)

- `width`:

  fig width

- `height`:

  fig height

#### Returns

path the file was written to (invisible)

------------------------------------------------------------------------

### `LFQDataPlotter$write_pdf()`

write figure to pdf

#### Usage

    LFQDataPlotter$write_pdf(fig, path_qc, fig_name, width = 7, height = 7)

#### Arguments

- `fig`:

  ggplot or ComplexHeatmap::Heatmap

- `path_qc`:

  path to write to

- `fig_name`:

  name of figure; `.pdf` is appended unless it ends in `.pdf`

- `width`:

  figure width

- `height`:

  figure height

#### Returns

path the file was written to

------------------------------------------------------------------------

### `LFQDataPlotter$clone()`

The objects of this class are cloneable with this method.

#### Usage

    LFQDataPlotter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> `protein_Id`: columns `proteinID`, `idtype2` united with separator '~'.
#> creating sampleName from file_name column
#> setup_analysis: dropped 6 unconfigured columns: average_prot_abundance, sd, mean, N, avg_peptide_abd, Replicate.
#> setup_analysis: added 36 rows with NA response for feature and sample combinations absent from the input.

lfqdata <- LFQData$new(
 istar$data,
 istar$config)
lfqplotter <- lfqdata$get_Plotter()

stopifnot(methods::is(lfqplotter$heatmap(), "Heatmap"))
stopifnot(methods::is(lfqplotter$heatmap_cor(), "Heatmap"))
stopifnot("ggplot" %in% class(lfqplotter$pca()))
#> PCA: removed 16 of 28 features with missing values. To keep all features, impute missing values first, e.g. AggregateLimpa$new(lfqdata, impute_only = TRUE)$aggregate().
stopifnot("plotly" %in%  class(lfqplotter$pca_plotly()))
#> PCA: removed 16 of 28 features with missing values. To keep all features, impute missing values first, e.g. AggregateLimpa$new(lfqdata, impute_only = TRUE)$aggregate().
tmp <- lfqplotter$boxplots()
stopifnot("ggplot" %in%  class(tmp$boxplot[[1]]))
stopifnot("ggplot" %in% class(lfqplotter$missingness_histogram()))

stopifnot(methods::is(lfqplotter$na_heatmap(), "Heatmap"))
class(lfqplotter$intensity_distribution_density())
#> [1] "ggplot2::ggplot" "ggplot"          "ggplot2::gg"     "S7_object"      
#> [5] "gg"             
class(lfqplotter$intensity_distribution_violin())
#> [1] "ggplot2::ggplot" "ggplot"          "ggplot2::gg"     "S7_object"      
#> [5] "gg"             
stopifnot(is.null(lfqplotter$pairs_smooth()))

stopifnot(class(lfqplotter$sample_correlation()) == "list")

stopifnot(methods::is(lfqplotter$raster(), "Heatmap"))
stopifnot("upset" == class(lfqplotter$upset_missing()))
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue at <https://github.com/hms-dbmi/UpSetR/issues>.
```
