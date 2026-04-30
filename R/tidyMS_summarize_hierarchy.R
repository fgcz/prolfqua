# Functions - summarize factors ----

#' Table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param file_name character — file name column
#' @param sample_name character — sample name column
#' @param factor_keys character vector — factor column names
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' xx <- table_factors(lfq$data_long(), lfq$file_name(), lfq$sample_name(), lfq$factor_keys())
#' xt <- xx |> dplyr::group_by(!!!rlang::syms(lfq$factor_keys())) |>
#'  dplyr::summarize(n = dplyr::n())
#' stopifnot(all(xt$n == 4))
#'
table_factors <- function(pdata, file_name, sample_name, factor_keys) {
  factors_table <- pdata |>
    dplyr::select(c(file_name, sample_name, factor_keys)) |>
    dplyr::distinct() |>
    arrange(!!sym(sample_name))
  return(factors_table)
}

#' Table of distinct factors with group sizes
#' @param pdata data.frame
#' @param file_name character — file name column
#' @param sample_name character — sample name column
#' @param factor_keys character vector — all factor column names
#' @param factor_keys_depth character vector — factor columns at current depth
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(istar$data, istar$config)
#' xx <- table_factors_size(lfq$data_long(), lfq$file_name(),
#'   lfq$sample_name(), lfq$factor_keys(), lfq$relevant_factor_keys())
#' stopifnot(all(xx$n == 4))
#'
table_factors_size <- function(pdata, file_name, sample_name, factor_keys, factor_keys_depth) {
  xx <- table_factors(pdata, file_name, sample_name, factor_keys)
  xx <- xx |>
    dplyr::group_by(dplyr::across(factor_keys_depth)) |>
    dplyr::summarize(n = dplyr::n())
  return(xx)
}

# Functions - summarize hierarchies

#' Count distinct elements for each level of hierarchy and isotope
#'
#' E.g. number of proteins, peptides, precursors in the dataset
#'
#' @param pdata data.frame
#' @param hierarchy_keys character vector — all hierarchy column names
#' @param isotope_label character — isotope label column name
#' @export
#' @keywords internal
#' @family summary
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb$data, bb$config)
#'
#' x <- hierarchy_counts(lfq$data_long(), lfq$hierarchy_keys(), lfq$isotope_label())
#' x$protein_Id
#' stopifnot(ncol(x) == length(lfq$hierarchy_keys()) + 1)
#' # select non existing protein
#' data0 <- lfq$data_long() |> dplyr::filter(protein_Id == "XYZ")
#' tmp <- hierarchy_counts(data0, lfq$hierarchy_keys(), lfq$isotope_label())
#' stopifnot(nrow(tmp) == 0)
hierarchy_counts <- function(pdata, hierarchy_keys, isotope_label = "isotopeLabel") {
  res <- pdata |>
    dplyr::group_by(across(all_of(isotope_label))) |>
    dplyr::summarise(across(all_of(hierarchy_keys), n_distinct))
  return(res)
}

#' Count distinct elements for each level of hierarchy per sample
#'
#' Provides wide, long, and plot views of hierarchy counts.
#'
#' @return An R6 class generator.
#' @export
#' @family summary
#' @examples
#' bb <- sim_lfq_data_protein_config()
#' counts <- HierarchyCountsSample$new(bb$data, bb$config)
#' head(counts$wide())
HierarchyCountsSample <- R6::R6Class(
  "HierarchyCountsSample",
  public = list(
    #' @field .summary summarised data frame
    .summary = NULL,
    #' @field .configuration AnalysisConfiguration object
    .configuration = NULL,
    #' @description Create a new HierarchyCountsSample
    #' @param pdata data frame
    #' @param configuration AnalysisConfiguration
    #' @param nr_children minimum number of children to include
    initialize = function(pdata, configuration, nr_children = 1) {
      self$.configuration <- configuration
      hierarchy <- configuration$hierarchy_keys()
      self$.summary <- pdata |>
        dplyr::filter(
          !is.na(!!rlang::sym(configuration$get_response())),
          !!rlang::sym(configuration$nr_children) >= .env$nr_children
        ) |>
        dplyr::group_by(across(all_of(c(configuration$isotope_label, configuration$sample_name)))) |>
        dplyr::summarise(across(all_of(hierarchy), dplyr::n_distinct))
    },
    #' @description Return wide-format summary
    wide = function() {
      self$.summary
    },
    #' @description Return long-format summary
    long = function() {
      self$.summary |>
        tidyr::pivot_longer(
          cols = -dplyr::all_of(c(self$.configuration$isotope_label, self$.configuration$sample_name)),
          names_to = "key",
          values_to = "nr"
        )
    },
    #' @description Return barplot of hierarchy counts
    plot = function() {
      long <- self$long()
      if (nrow(long) == 0) {
        return(NULL)
      }
      nudgeval <- -mean(long$nr) * 0.05
      ggplot2::ggplot(long, ggplot2::aes(x = !!rlang::sym(self$.configuration$sample_name), y = .data$nr)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "white") +
        ggplot2::facet_wrap(~key, scales = "free_y", ncol = 1) +
        ggplot2::geom_text(ggplot2::aes(label = .data$nr), nudge_y = nudgeval, angle = 65) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    }
  )
)

#' Hierarchy counts per sample
#'
#' @param pdata data.frame
#' @param configuration AnalysisConfiguration
#' @param nr_children minimum number of children
#' @return \code{\link{HierarchyCountsSample}} R6 object
#' @export
#' @family summary
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#'
#' config <- bb$config
#' data <- bb$data
#' res <- hierarchy_counts_sample(data, config, nr_children = 1)
#' x <- res$long()
#' # filters on peptide level
#' res <- hierarchy_counts_sample(data, config, nr_children = 2)
#' x2 <- res$long()
#' # filters on protein level based on peptide count
#' bb <- prolfqua::sim_lfq_data_protein_config()
#' res <- hierarchy_counts_sample(bb$data, bb$config, nr_children = 2)
#' x1 <- res$wide()
#' res <- hierarchy_counts_sample(bb$data, bb$config, nr_children = 1)
#' x2 <- res$wide()
#' x1$nr_children <- 2
#' x2$nr_children <- 1
#' xl <- dplyr::bind_rows(x1, x2)
#'
#' xl$nr_children |> table()
#' nudgeval <-  -mean(xl$protein_Id) * 0.05
#' ggplot2::ggplot(xl,
#'   ggplot2::aes(x = sampleName, y = protein_Id, fill = as.character(nr_children))) +
#'  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge())
#'
hierarchy_counts_sample <- function(
  pdata,
  configuration,
  nr_children = 1
) {
  HierarchyCountsSample$new(pdata, configuration, nr_children = nr_children)
}


#' Summarize hierarchy counts
#'
#' E.g compute number of peptides for each protein
#'
#' @export
#' @keywords internal
#' @family summary
#' @param pdata data.frame
#' @param hierarchy_keys character vector — all hierarchy column names
#' @param isotope_label character — isotope label column name
#' @param hierarchy character vector — hierarchy level to group by (default hierarchy_keys)
#' @param factors character vector — factor columns to include
#'
#' @examples
#'
#' bb <- sim_lfq_data_peptide_config()
#' lfq <- LFQData$new(bb$data, bb$config)
#' summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), lfq$isotope_label())
#' summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), lfq$isotope_label(),
#'  factors = character())
#' summarize_hierarchy(lfq$data_long(), lfq$hierarchy_keys(), lfq$isotope_label(),
#'  hierarchy = lfq$relevant_hierarchy_keys())
#'
summarize_hierarchy <- function(
  pdata,
  hierarchy_keys,
  isotope_label = "isotopeLabel",
  hierarchy = hierarchy_keys,
  factors = character()
) {
  all_hierarchy <- c(isotope_label, hierarchy_keys)

  precursor <- pdata |> dplyr::select(dplyr::all_of(c(factors, all_hierarchy))) |> dplyr::distinct()
  x3 <- precursor |>
    dplyr::group_by(across(all_of(c(factors, hierarchy)))) |>
    dplyr::summarize(across(all_of(base::setdiff(all_hierarchy, hierarchy)), list(n = dplyr::n_distinct)))
  return(x3)
}
