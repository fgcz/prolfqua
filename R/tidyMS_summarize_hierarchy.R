# Functions - summarize factors ----

#' Table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param configuration AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#'
#'
#' xx <- table_factors(istar$data,istar$config )
#' xx
#' xt <- xx |> dplyr::group_by(!!!rlang::syms(istar$config$factor_keys())) |>
#'  dplyr::summarize(n = dplyr::n())
#' xt
#' stopifnot(all(xt$n == 4))
#'
table_factors <- function(pdata, configuration) {
  factorsTab <- pdata |>
    dplyr::select(c(configuration$fileName, configuration$sampleName, configuration$factor_keys())) |>
    dplyr::distinct() |>
    arrange(!!sym(configuration$sampleName))
  return(factorsTab)
}

#' Table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param configuration AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#'
#' istar <- sim_lfq_data_peptide_config()
#'
#'
#' xx <- table_factors_size(istar$data,istar$config )
#' stopifnot(all(xx$n == 4))
#'
table_factors_size <- function(pdata, configuration) {
  xx <- table_factors(pdata, configuration)
  xx <- xx |> dplyr::group_by(dplyr::across(configuration$factor_keys_depth())) |> dplyr::summarize(n = dplyr::n())
  return(xx)
}

# Functions - summarize hierarchies

#' Count distinct elements for each level of hierarchy and istope
#'
#' E.g. number of proteins, peptides, precursors in the dataset
#'
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @export
#' @keywords internal
#' @family summary
#' @examples
#'
#' bb <- prolfqua::sim_lfq_data_peptide_config()
#'
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' x <- hierarchy_counts(data, config)
#' x$protein_Id
#' stopifnot(ncol(x) == length(config$hierarchy_keys()) + 1)
#' # select non existing protein
#' data0 <- data |> dplyr::filter( protein_Id == "XYZ")
#' tmp <- hierarchy_counts(data0, config)
#' stopifnot(nrow(tmp) == 0)
hierarchy_counts <- function(pdata, config) {
  hierarchy <- config$hierarchy_keys()
  res <- pdata |>
    dplyr::group_by(across(all_of(config$isotopeLabel))) |>
    dplyr::summarise(across(all_of(hierarchy), n_distinct))

  return(res)
}

#' Count distinct elements for each level of hierarchy per sample
#'
#' @export
#' @param pdata data.frame
#' @param configuration \code{\link{AnalysisConfiguration}}
#' @param nr_children integer, minimum number of children required
#'
#' @keywords internal
#' R6 class for hierarchy counts per sample
#'
#' Provides wide, long, and plot views of hierarchy counts.
#'
#' @export
#' @family summary
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
        dplyr::group_by(across(all_of(c(configuration$isotopeLabel, configuration$sampleName)))) |>
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
          cols = -dplyr::all_of(c(self$.configuration$isotopeLabel, self$.configuration$sampleName)),
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
      ggplot2::ggplot(long, ggplot2::aes(x = !!rlang::sym(self$.configuration$sampleName), y = .data$nr)) +
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
#' @param config AnalysisConfiguration
#' @param hierarchy for which hierarchy level (default up to hierarchy depth)
#' @param factors which factors to include
#'
#' @examples
#'
#'
#'
#' bb <- sim_lfq_data_peptide_config()
#' data <- bb$data
#' configur <- bb$config
#' summarize_hierarchy(data, configur)
#' summarize_hierarchy(data, configur, factors = character())
#'
#' summarize_hierarchy(data, configur,
#'  hierarchy = configur$hierarchy_keys_depth() )
#' summarize_hierarchy(data, configur,
#'  hierarchy = NULL, factors = configur$factor_keys_depth() )
#' configur$hierarchyDepth = 1
#' summarize_hierarchy(data, configur,
#'  factors = configur$factor_keys_depth())
#' configur$hierarchyDepth = 2
#' summarize_hierarchy(data, configur)
#' configur$hierarchyDepth = 3
#' summarize_hierarchy(data, configur )
#' configur$hierarchyDepth = 4
#' summarize_hierarchy(data, configur )
#'
summarize_hierarchy <- function(pdata, config, hierarchy = config$hierarchy_keys_depth(), factors = character()) {
  all_hierarchy <- c(config$isotopeLabel, config$hierarchy_keys())

  precursor <- pdata |> dplyr::select(dplyr::all_of(c(factors, all_hierarchy))) |> dplyr::distinct()
  x3 <- precursor |>
    dplyr::group_by(across(all_of(c(factors, hierarchy)))) |>
    dplyr::summarize(across(all_of(base::setdiff(all_hierarchy, hierarchy)), list(n = dplyr::n_distinct)))
  return(x3)
}
