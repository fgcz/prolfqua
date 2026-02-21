# LFQDataPlotter —- Create various visualization of the LFQdata

LFQDataPlotter —- Create various visualization of the LFQdata

LFQDataPlotter —- Create various visualization of the LFQdata

## See also

[`plot_heatmap_cor`](https://wolski.github.io/prolfqua/reference/plot_heatmap_cor.md)

[`plot_pca`](https://wolski.github.io/prolfqua/reference/plot_pca.md)

Other LFQData:
[`LFQData`](https://wolski.github.io/prolfqua/reference/LFQData.md),
[`LFQDataAggregator`](https://wolski.github.io/prolfqua/reference/LFQDataAggregator.md),
[`LFQDataImp`](https://wolski.github.io/prolfqua/reference/LFQDataImp.md),
[`LFQDataStats`](https://wolski.github.io/prolfqua/reference/LFQDataStats.md),
[`LFQDataSummariser`](https://wolski.github.io/prolfqua/reference/LFQDataSummariser.md),
[`LFQDataToSummarizedExperiment()`](https://wolski.github.io/prolfqua/reference/LFQDataToSummarizedExperiment.md)

## Public fields

- `lfq`:

  LFQData object

- `prefix`:

  prefix to figure names when writing, e.g. protein\_

- `file_paths_pdf`:

  with paths to figures

- `file_paths_html`:

  with paths to figures

## Methods

### Public methods

- [`LFQDataPlotter$new()`](#method-LFQDataPlotter-new)

- [`LFQDataPlotter$raster()`](#method-LFQDataPlotter-raster)

- [`LFQDataPlotter$heatmap()`](#method-LFQDataPlotter-heatmap)

- [`LFQDataPlotter$heatmap_cor()`](#method-LFQDataPlotter-heatmap_cor)

- [`LFQDataPlotter$pca()`](#method-LFQDataPlotter-pca)

- [`LFQDataPlotter$pca_plotly()`](#method-LFQDataPlotter-pca_plotly)

- [`LFQDataPlotter$boxplots()`](#method-LFQDataPlotter-boxplots)

- [`LFQDataPlotter$missigness_histogram()`](#method-LFQDataPlotter-missigness_histogram)

- [`LFQDataPlotter$NA_heatmap()`](#method-LFQDataPlotter-NA_heatmap)

- [`LFQDataPlotter$intensity_distribution_density()`](#method-LFQDataPlotter-intensity_distribution_density)

- [`LFQDataPlotter$intensity_distribution_violin()`](#method-LFQDataPlotter-intensity_distribution_violin)

- [`LFQDataPlotter$pairs_smooth()`](#method-LFQDataPlotter-pairs_smooth)

- [`LFQDataPlotter$sample_correlation()`](#method-LFQDataPlotter-sample_correlation)

- [`LFQDataPlotter$upset_missing()`](#method-LFQDataPlotter-upset_missing)

- [`LFQDataPlotter$write_boxplots()`](#method-LFQDataPlotter-write_boxplots)

- [`LFQDataPlotter$write_pltly()`](#method-LFQDataPlotter-write_pltly)

- [`LFQDataPlotter$write_pdf()`](#method-LFQDataPlotter-write_pdf)

- [`LFQDataPlotter$write()`](#method-LFQDataPlotter-write)

- [`LFQDataPlotter$clone()`](#method-LFQDataPlotter-clone)

------------------------------------------------------------------------

### Method `new()`

create LFQDataPlotter

#### Usage

    LFQDataPlotter$new(lfqdata, prefix = "ms_")

#### Arguments

- `lfqdata`:

  LFQData

- `prefix`:

  will be prepended to outputs written

------------------------------------------------------------------------

### Method `raster()`

plot intensities in raster

#### Usage

    LFQDataPlotter$raster(
      arrange = c("mean", "var"),
      not_na = FALSE,
      rownames = FALSE
    )

#### Arguments

- `arrange`:

  arrange by either mean or var

- `not_na`:

  TRUE arrange by number of NA's, FALSE by arrange by intensity

- `rownames`:

  show rownames (default FALSE - do not show.)

#### Returns

ggplot

------------------------------------------------------------------------

### Method [`heatmap()`](https://rdrr.io/r/stats/heatmap.html)

heatmap of intensities - columns are samples, rows are proteins or
peptides.

The abundances of each protein (row) are z-scored. Afterward, the mean
abundance for each protein is zero, and the standard variation is one.
z-scoring allows to compare (cluster) the proteins according to the
difference in the expression in the samples. Without the z-scoring, the
proteins would group according to their abundance, e.g., high abundant
proteins would be one cluster.

#### Usage

    LFQDataPlotter$heatmap(na_fraction = 0.3, rownames = FALSE)

#### Arguments

- `na_fraction`:

  max fraction of NA's per row

- `rownames`:

  show rownames (default FALSE - do not show.)

#### Returns

pheatmap

------------------------------------------------------------------------

### Method `heatmap_cor()`

heatmap of sample correlations.

The Spearman correlation among all samples is computed. Then the
euclidean distance is used to compute the distances.

#### Usage

    LFQDataPlotter$heatmap_cor()

#### Returns

pheatmap

------------------------------------------------------------------------

### Method `pca()`

PCA plot

A PCA is applied and the first and second principal component are shown.
Features with missing values are removed. For PCA with all features,
impute first using
[`LFQDataImp`](https://wolski.github.io/prolfqua/reference/LFQDataImp.md).

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

### Method `pca_plotly()`

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

### Method `boxplots()`

boxplots for all proteins

#### Usage

    LFQDataPlotter$boxplots(facet = TRUE)

#### Arguments

- `facet`:

  enable facet wrap if hierarchy_depth less then hierarchy lenght.

#### Returns

tibble with column boxplots containing ggplot objects

------------------------------------------------------------------------

### Method [`missigness_histogram()`](https://wolski.github.io/prolfqua/reference/missigness_histogram.md)

histogram of intensities given number of missing in conditions

#### Usage

    LFQDataPlotter$missigness_histogram()

#### Returns

ggplot

------------------------------------------------------------------------

### Method `NA_heatmap()`

heatmap of features with missing values

#### Usage

    LFQDataPlotter$NA_heatmap()

#### Returns

ggplot

------------------------------------------------------------------------

### Method `intensity_distribution_density()`

density distribution of intensities

#### Usage

    LFQDataPlotter$intensity_distribution_density(legend = TRUE)

#### Arguments

- `legend`:

  show legend TRUE, FALSE do not show.

#### Returns

ggplot

------------------------------------------------------------------------

### Method `intensity_distribution_violin()`

Violinplot showing distribution of intensities in all samples

#### Usage

    LFQDataPlotter$intensity_distribution_violin()

#### Returns

ggplot

------------------------------------------------------------------------

### Method [`pairs_smooth()`](https://wolski.github.io/prolfqua/reference/pairs_smooth.md)

pairsplot of intensities

#### Usage

    LFQDataPlotter$pairs_smooth(max = 10)

#### Arguments

- `max`:

  maximal number of samples to show

#### Returns

NULL

------------------------------------------------------------------------

### Method `sample_correlation()`

plot of sample correlations

#### Usage

    LFQDataPlotter$sample_correlation()

#### Returns

NULL

------------------------------------------------------------------------

### Method `upset_missing()`

upset plot based on presence absence information

#### Usage

    LFQDataPlotter$upset_missing()

#### Returns

plot

------------------------------------------------------------------------

### Method `write_boxplots()`

write boxplots to file

#### Usage

    LFQDataPlotter$write_boxplots(path_qc, filename = NULL, width = 6, height = 6)

#### Arguments

- `path_qc`:

  path to write to

- `filename`:

  file to write into

- `width`:

  fig width

- `height`:

  fig height

------------------------------------------------------------------------

### Method `write_pltly()`

write pltly figures to path_qc

#### Usage

    LFQDataPlotter$write_pltly(fig, path_qc, fig_name)

#### Arguments

- `fig`:

  pltly figure

- `path_qc`:

  path to write to

- `fig_name`:

  file name (without extension)

#### Returns

path the file was written to.

------------------------------------------------------------------------

### Method `write_pdf()`

write figure to pdf

#### Usage

    LFQDataPlotter$write_pdf(fig, path_qc, fig_name, width = 7, height = 7)

#### Arguments

- `fig`:

  ggplot or pheatmap

- `path_qc`:

  path to write to

- `fig_name`:

  name of figure (no extension)

- `width`:

  figure width

- `height`:

  figure height

#### Returns

path the file was written to

------------------------------------------------------------------------

### Method [`write()`](https://rdrr.io/r/base/write.html)

write heatmaps and pca plots to files

#### Usage

    LFQDataPlotter$write(path_qc)

#### Arguments

- `path_qc`:

  path to write to

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    LFQDataPlotter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- sim_lfq_data_peptide_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done

lfqdata <- LFQData$new(
 istar$data,
 istar$config)
lfqplotter <- lfqdata$get_Plotter()

stopifnot(class(lfqplotter$heatmap()) == "pheatmap")
stopifnot(class(lfqplotter$heatmap_cor()) == "pheatmap")
stopifnot("ggplot" %in% class(lfqplotter$pca()))
#> PCA: removed 16 of 28 features with missing values. For PCA with all features, impute first using LFQDataImp.
#> Joining with `by = join_by(sampleName)`
stopifnot("plotly" %in%  class(lfqplotter$pca_plotly()))
#> PCA: removed 16 of 28 features with missing values. For PCA with all features, impute first using LFQDataImp.
#> Joining with `by = join_by(sampleName)`
tmp <- lfqplotter$boxplots()
stopifnot("ggplot" %in%  class(tmp$boxplot[[1]]))
stopifnot("ggplot" %in% class(lfqplotter$missigness_histogram()))
#> Warning: >>>> deprecated! <<<< 
#> 
#>           use summarize_stats_factors instead.
#> completing cases
#> isotopeLabel ~ group_

stopifnot(class(lfqplotter$NA_heatmap()) == "pheatmap")
#> rows with NA's: 16; all rows :28
class(lfqplotter$intensity_distribution_density())
#> [1] "ggplot2::ggplot" "ggplot"          "ggplot2::gg"     "S7_object"      
#> [5] "gg"             
class(lfqplotter$intensity_distribution_violin())
#> [1] "ggplot2::ggplot" "ggplot"          "ggplot2::gg"     "S7_object"      
#> [5] "gg"             
stopifnot(is.null(lfqplotter$pairs_smooth()))

stopifnot(class(lfqplotter$sample_correlation()) == "list")

stopifnot(class(lfqplotter$raster()) == "pheatmap")
stopifnot("upset" == class(lfqplotter$upset_missing()))
#> completing cases
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.
stopifnot(class(prolfqua::plot_sample_correlation(istar$data, istar$config)) == "list")

```
