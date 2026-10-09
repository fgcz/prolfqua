#' @importFrom forcats fct_relevel
#' @importFrom dplyr across all_of anti_join arrange bind_cols bind_rows case_when count desc distinct filter group_by
#' @importFrom dplyr inner_join left_join mutate nest_by
#' @importFrom dplyr rename select starts_with summarize ungroup
#'
#' @importFrom ggplot2 aes element_text facet_grid ggplot ggtitle geom_boxplot
#' @importFrom ggplot2 geom_line geom_violin guides
#' @importFrom ggplot2 stat_summary scale_x_continuous scale_y_continuous theme
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom graphics par rect hist
#' @importFrom graphics abline pairs plot text
#' @importFrom grDevices colorRampPalette dev.off pdf png rainbow
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom plotly ggplotly
#' @importFrom purrr map map2 map2_dbl map_lgl map_chr map_dbl reduce map_if map_dfc map_int map_df
#' @importFrom rlang := sym syms .data
#' @importFrom stats as.formula cor
#' @importFrom stats .lm.fit glm model.matrix residuals
#' @importFrom tidyr nest nesting separate_rows unite unnest separate
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble add_column as_tibble column_to_rownames tibble
#' @importFrom stats anova coef coefficients confint cor.test delete.response
#'   df.residual ecdf formula lm mad median medpolish
#' @importFrom stats na.omit p.adjust pbeta power.t.test prcomp predict pt qt
#'   quantile sd setNames sigma terms update vcov
#' @importFrom stats approx df fisher.test ftable lm.fit pchisq pf pnorm qf rgeom rlnorm rnorm uniroot
#' @importFrom utils combn data head read.csv tail unzip
#' @importFrom limma lmFit contrasts.fit eBayes topTable
NULL


# Suppress R CMD check NOTEs for variables used in NSE (dplyr/ggplot2)
utils::globalVariables(c(
  ".",
  ".env",
  "Freq",
  "Group",
  "PC",
  "Replicate",
  "group",
  "isSingular",
  "meanAbundance",
  "moderated.var.post",
  "nrPeptides",
  "reference_mean",
  "reference_median",
  "rhs",
  "subject",
  "y"
))

#' Internal Functions by category
#' @family aggregation
#' @family configuration
#' @family modelling
#' @family plotting
#' @family preprocessing
#' @family stats
#' @family summary
#' @family utilities
#' @name INTERNAL_FUNCTIONS_BY_FAMILY
#'
#' @return A documentation topic.
NULL

#' Package Data
#' @family data
#' @name PACKAGE_DATA
#' @return A documentation topic.
NULL


#' Factor annotation of 36 samples
#'
#' Sample annotation table with two factors (`factor_A`: MI / MINOCA,
#' `factor_B`: sampling time T0 / T150 / T300), used to demonstrate how
#' contrasts are derived from the factor levels. See `inst/script/README.md`
#' for its origin.
#'
#' @format A tibble with 36 rows and the columns `Relative Path`, `Name`,
#'   `Group`, `factor_A`, `factor_B`.
#' @family data
#' @docType data
#' @keywords internal
x5463yzwer453bbb <- NULL

#' Reference input for the robust summarisation check (87 percent)
#'
#' log2 peptide areas of one protein, used by the tests of the robust
#' (medpolish/rlm) protein summarisation. See `inst/script/README.md`.
#'
#' @format A data frame with 20 rows and the columns `log2Area`, `peptide_Id`,
#'   `sampleName`.
#' @family data
#' @docType data
#' @keywords internal
data_checksummarizationrobust87 <- NULL

#' Reference input for the robust summarisation check (69 percent)
#'
#' log2 peptide areas of one protein, used by the tests of the robust
#' (medpolish/rlm) protein summarisation. See `inst/script/README.md`.
#'
#' @format A data frame with 20 rows and the columns `log2Area`, `peptide_Id`,
#'   `sampleName`.
#' @family data
#' @docType data
#' @keywords internal
data_checksummarizerobust69 <- NULL

#' IonStar peptide-level example data
#'
#' Peptide intensities of 200 proteins from the IonStar human cell lysate with
#' E. coli spike-in dilution series (dilutions a to e, 4 replicates each),
#' processed with MaxQuant. An R6 object holding the data, its
#' \code{\link{AnalysisConfiguration}}, and the methods \code{filtered()},
#' \code{normalized()} and \code{subset_normalized()}. See
#' `inst/script/README.md` for its origin.
#'
#' @format An R6 object of class `IonstarData`; `$data` is a data frame with
#'   26520 rows.
#' @source PRIDE PXD003881
#' @family data
#' @docType data
#' @keywords internal
data_ionstar <- NULL
