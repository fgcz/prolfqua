# AnalysisConfiguration ----
#'
#' Analysis Configuration
#' @description
#' Analysis Configuration
#'
#' @family configuration
#' @export
#'
AnalysisConfiguration <- R6::R6Class(
  "AnalysisConfiguration",
  public = list(
    #' @field sep separator to use when uniting columns is necessary
    sep = "~",
    #' @field table AnalysisTableAnnotation
    table = NULL,
    #' @field parameter AnalysisParameter
    parameter = NULL,
    #' @description
    #' create AnalysisConfiguration
    #' @param analysisTableAnnotation instance of AnalysisTableAnnotation
    #' @param analysisParameter instance of AnalysisParameter
    initialize = function(analysisTableAnnotation, analysisParameter = AnalysisParameters$new()){
      self$table <- analysisTableAnnotation
      self$parameter <- analysisParameter
    }
  )
)

#' Make reduced hierarchy configuration
#' @export
#' @keywords internal
#' @param config AnalysisConfiguration
#' @param workIntensity work intensity column
#' @param hierarchy new reduced hierarchy
#' @family configuration
#' @return AnalysisConfiguration with reduced hieararchy
#' @examples
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' make_reduced_hierarchy_config(conf,
#'  "testintensity",
#'  conf$table$hierarchy[1:2])
#'
make_reduced_hierarchy_config <- function(config, workIntensity , hierarchy ){
  newConfig <- config$clone(deep = TRUE)
  newConfig$table$hierarchy = hierarchy
  newConfig$table$workIntensity = workIntensity
  return(newConfig)
}

#' create interaction column from factors
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#' xx <- data.frame(A = c("a","a","a"), B = c("d","d","e"))
#' # debug(make_interaction_column)
#' x <- make_interaction_column(xx, c("B","A"))
#' x <- make_interaction_column(xx, c("A"))
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' config$table$factor_keys()
#' config$table$factorDepth <- 1
#' make_interaction_column(analysis,
#'    config$table$factor_keys_depth())
#'
make_interaction_column <- function(data, columns, sep="."){
  intr <- dplyr::select(data, columns)
  intr <- purrr::map_dfc(intr, factor)
  names(columns) <- columns
  newlev <- purrr::map2(columns, intr, function(x,y){paste0(x,levels(y))})
  intr <- purrr::map2_dfc(columns, intr, paste0)
  intr <- purrr::map2_dfc(intr , newlev, forcats::fct_relevel)

  colnames(intr) <- paste0("interaction_",columns)
  colname <- "interaction"
  data <- data |> dplyr::mutate(!!colname := interaction(intr, sep = sep))
  return(data)
}


# Functions - Configuration ----
#' Extract all value slots in an R6 object
#' @param r6class r6 class
#' @keywords internal
#' @family configuration
#' @export
R6_extract_values <- function(r6class){
  tmp <- sapply(r6class, class)
  slots <- tmp[ !tmp %in% c("environment", "function")]
  res <- list()
  for (i in names(slots)) {
    if ("R6" %in% class(r6class[[i]])) {
      res[[i]]  <- R6_extract_values(r6class[[i]])
    }else{
      res[[i]] <- r6class[[i]]
    }
  }
  return(res)
}


#' Setup a tidy table compatible with a \code{\link{AnalysisConfiguration}}
#'
#' Extracts columns relevant for a configuration from a data frame
#' and create new columns e.g. sampleName column etc.
#' @param data data.frame
#' @param config AnlalysisConfiguration
#' @param cc complete cases default TRUE
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#' skylineconfig <- create_config_Skyline(isotopeLabel = "Isotope.Label.Type",
#'  ident_qValue = "Detection.Q.Value")
#' skylineconfig$table$factors[["Time"]] = "Sampling.Time.Point"
#'
#' sample_analysis <- setup_analysis(prolfqua_data('data_skylinePRMSample_A')$data, skylineconfig)
#'
setup_analysis <- function(data, configuration, cc = TRUE,  from_factors = FALSE){
  configuration <- configuration$clone(deep = TRUE)
  table <- configuration$table

  # extract hierarchy columns
  for (i in seq_along(table$hierarchy))
  {
    data <- tidyr::unite(data, UQ(sym(table$hierarchy_keys()[i])), table$hierarchy[[i]],remove = FALSE, sep = configuration$sep)
  }
  data <- dplyr::select(data , -dplyr::one_of(dplyr::setdiff(unlist(table$hierarchy), table$hierarchy_keys() )))

  # extract factors
  if ( length(table$factors) == 0) {
    stop("No factors (explanatory variables) specified in the AnalysisTableConfiguration.\n",
         'Pleases use table$factors["Condition"] = "columnName".\n',
         'where Condition is the new name of the variable and\n',
         'columnName is the name of the column containing the varible.')
  }
  for (i in seq_along(table$factors))
  {
    if ( length(table$factors[[i]]) > 1) {
      data <- tidyr::unite(data, UQ(sym(table$factor_keys()[i])), table$factors[[i]],remove = FALSE, sep = configuration$sep)
    } else {
      newname <- table$factor_keys()[i]
      data <- dplyr::mutate(data, !!newname := as.character(!!sym(table$factors[[i]])))
    }
  }

  sampleName <- table$sampleName

  if (from_factors & !sampleName  %in% names(data)) {
    message("creating sampleName from factor columns")
    data <- data |>  tidyr::unite(
      UQ(sym(sampleName)) ,
      unique(unlist(table$factors)), remove = TRUE , sep = configuration$sep) |>
      dplyr::select(sampleName, table$fileName) |> dplyr::distinct() |>
      dplyr::mutate_at(sampleName, function(x){ x <- make.unique( x, sep = configuration$sep )}) |>
      dplyr::inner_join(data, by = table$fileName)
  } else if (!sampleName  %in% names(data)) {
    message("creating sampleName from fileName column")
    data[[table$sampleName]] <- tools::file_path_sans_ext( basename(data[[table$fileName]]) )
  }else {
    message(sampleName, " already exists")
  }

  data <- data |> dplyr::select(-dplyr::one_of(dplyr::setdiff(unlist(table$factors), table$factor_keys())))

  # Make implicit NA's explicit
  if (!(table$isotopeLabel %in% colnames(data))) {
    warning("no isotopeLabel column specified in the data, adding column automatically and setting to 'light'.")
    data[[table$isotopeLabel]] <- "light"
  }
  if (!(table$ident_qValue %in%  colnames(data))) {
    warning("no qValue column specified in the data. Creating column and setting qValues to 0.")
    data[[table$ident_qValue]] <- 0
  }

  # TODO add better warning....
  data <- data |> dplyr::select(c(table$id_vars(),table$value_vars()))

  txd <- data |> group_by(!!!syms(c(table$fileName, table$hierarchy_keys(), table$isotopeLabel))) |>
    summarize(n = n())
  if (length(table(txd$n)) > 1) {
    str <- paste("There is more than ONE observations for each : ", paste( table$hierarchy_keys(), collapse = ", "), ",\n",
                 "and sample : ", table$sampleName, "; (filename) : ", table$fileName, "\n")
    warning(str)
    warning("Please inspect the returned dataframe. Check for rows where n > 1\n e.g. res |> dplyr::filter(n > 1)")
    return(txd)
  }
  #tmp <- prolfqua::tidy_to_wide_config(data, configuration)
  #message("nr rows and nr columns")
  #message(paste(dim(tmp$data),collapse = ", "))

  if (cc) {
    data <- complete_cases( data , configuration)
  }


  return( data )
}

#' Separates hierarchy columns into starting columns
#'
#'
#' @export
#' @seealso \code{\link{lfq_write_table}}
#' @param data data.frame
#' @param config AnlalysisConfiguration
#' @family configuration
#' @keywords internal
#' @examples
#'
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' dt <- separate_hierarchy(bb$data, bb$config)
#' setdiff(colnames(dt) ,colnames(bb$data))
#' stopifnot(ncol(dt) >= ncol(bb$data))
#'
separate_hierarchy <- function(data, config){
  for (hkey in config$table$hierarchy_keys_depth()) {
    if (length(config$table$hierarchy[[hkey]]) == 1 && hkey == config$table$hierarchy[[hkey]]) {

    }else {
      data <- data |> tidyr::separate( hkey, config$table$hierarchy[[hkey]], sep = config$sep, remove = FALSE)
    }
  }
  return(data)
}

#'
#' Separates factor columns into starting columns
#'
#' @param data data.frame
#' @param config AnlalysisConfiguration
#' @export
#'
#' @keywords internal
#' @family configuration
#' @examples
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' dt <- separate_factors(bb$data, bb$config)
#' setdiff(colnames(dt), colnames(bb$data))
#' stopifnot(ncol(bb$data) < ncol(dt))
#'
separate_factors <- function(data, config) {
  for (fkey in config$table$factor_keys()) {
    data <- data |> tidyr::separate( fkey, config$table$factors[[fkey]], sep = config$sep, remove = FALSE)
  }
  return(data)
}


#' Complete cases
#'
#' The tidy table does not need to contain missing data.
#' This function re-establishes the missing observations in a sample.
#'
#' @export
#' @param pdata data.frame
#' @param config AnlalysisConfiguration
#' @keywords internal
#' @family configuration
#' @examples
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#'
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' xx <- complete_cases(data, config)
#' stopifnot(nrow(data) <= nrow(xx))
#'
complete_cases <- function(pdata, config) {
  message("completing cases")
  fkeys <- c(config$table$fileName,config$table$sampleName, config$table$factor_keys())
  hkeys <- c(config$table$isotopeLabel, config$table$hierarchy_keys())
  res <- tidyr::complete(
    pdata,
    tidyr::nesting(!!!syms(fkeys)),
    tidyr::nesting(!!!syms(hkeys))
  )
  return(res)
}


#' Sample subset of proteins/peptides/precursors
#' @param n size of sample
#' @param pdata tidy table
#' @param config \code{\link{AnalysisConfiguration}}
#' @export
#' @keywords internal
#' @family configuration
#'
sample_subset <- function(size, pdata, config){
  hk <- config$table$hierarchy_keys_depth()
  message("Sampling ", size, paste(hk, collapse = "," ) )
  hkdf <- pdata |> select(all_of(hk)) |> distinct() |> sample_n(size = size)
  sdata <- inner_join(hkdf, pdata)
  return(sdata)
}

# Functions - summary ----

# Functions - summarize factors ----

#' Table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' xx <- table_factors(data,config )
#' xx
#' xt <- xx |> dplyr::group_by(!!!rlang::syms(config$table$factor_keys())) |>
#'  dplyr::summarize(n = dplyr::n())
#' stopifnot(all(xt$n == 1))
#'
table_factors <- function(pdata, configuration){
  factorsTab <- pdata |> dplyr::select(c(configuration$table$fileName, configuration$table$sampleName, configuration$table$factor_keys())) |>
    dplyr::distinct() |>
    arrange(!!sym(configuration$table$sampleName))
  return(factorsTab)
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
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' x <- hierarchy_counts(data, config)
#' x$protein_Id
#' stopifnot(ncol(x) == length(config$table$hierarchy_keys()) + 1)
#' # select non existing protein
#' data <- data |> dplyr::filter( protein_Id == "XYZ")
#' tmp <- hierarchy_counts(data, config)
#' stopifnot(nrow(tmp) == 0)
hierarchy_counts <- function(pdata, config){
  hierarchy <- config$table$hierarchy_keys()
  res <- pdata |>
    dplyr::group_by_at(config$table$isotopeLabel) |>
    dplyr::summarise_at( hierarchy, n_distinct )

  return(res)
}

#' Count distinct elements for each level of hierarchy per sample
#'
#' @export
#' @param pdata data.frame
#' @param config \code{\link{AnalysisConfiguration}}
#'
#' @keywords internal
#' @family summary
#' @examples
#'
#'
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' res <- hierarchy_counts_sample(data, config)
#' res()
#' res("long")
#' res("plot")
hierarchy_counts_sample <- function(pdata,
                                    configuration)
{
  hierarchy <- configuration$table$hierarchy_keys()
  summary <- pdata |> dplyr::filter(!is.na(!!sym(configuration$table$get_response() ))) |>
    dplyr::group_by_at(c(configuration$table$isotopeLabel, configuration$table$sampleName)) |>
    dplyr::summarise_at( hierarchy, n_distinct )

  res <- function(value = c("wide", "long", "plot")){
    value <- match.arg(value)
    if (value == "wide") {
      return(summary)
    }else{
      long <- summary |> tidyr::gather("key",
                                       "nr",-dplyr::one_of(configuration$table$isotopeLabel,
                                                           configuration$table$sampleName))
      if (value == "long") {
        return(long)
      }else if (value == "plot") {
        nudgeval <-  -mean(long$nr) * 0.05
        # TODO(WEW) check potential problem with sampleName
        ggplot2::ggplot(long, aes(x = !!sym(configuration$table$sampleName), y = .data$nr)) +
          ggplot2::geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "white") +
          ggplot2::facet_wrap( ~ key, scales = "free_y", ncol = 1) +
          ggplot2::geom_text(aes(label = .data$nr), nudge_y = nudgeval, angle = 65) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      }
    }
  }
  return(res)
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
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' summarize_hierarchy(data, configur)
#' summarize_hierarchy(data, configur, factors = character())
#'
#' summarize_hierarchy(data, configur,
#'  hierarchy = configur$table$hierarchy_keys_depth() )
#' summarize_hierarchy(data, configur,
#'  hierarchy = NULL, factors = configur$table$factor_keys_depth() )
#' configur$table$hierarchyDepth = 1
#' summarize_hierarchy(data, configur,
#'  factors = configur$table$factor_keys_depth())
#' configur$table$hierarchyDepth = 2
#' summarize_hierarchy(data, configur)
#' configur$table$hierarchyDepth = 3
#' summarize_hierarchy(data, configur )
#' configur$table$hierarchyDepth = 4
#' summarize_hierarchy(data, configur )
#'
summarize_hierarchy <- function(pdata,
                                config,
                                hierarchy = config$table$hierarchy_keys_depth(),
                                factors = character())
{
  all_hierarchy <- c(config$table$isotopeLabel, config$table$hierarchy_keys() )

  precursor <- pdata |> dplyr::select(factors, all_hierarchy) |> dplyr::distinct()
  x3 <- precursor |> dplyr::group_by_at(c(factors, hierarchy)) |>
    dplyr::summarize_at( setdiff(all_hierarchy, hierarchy),
                         list( n = dplyr::n_distinct))
  return(x3)
}


# Functions - Handling isotopes ----

#' Spreads isotope label heavy and light into two columns
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' x<-spread_response_by_IsotopeLabel(data,configur)
#'
#' bb <- prolfqua_data('data_skylineSRM_HL_A')
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' bb <- prolfqua_data('data_skylineSRM_HL_A')
#' conf <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' x <- spread_response_by_IsotopeLabel(analysis, conf)
#'
spread_response_by_IsotopeLabel <- function(resData, config){
  table <- config$table
  id_vars <- table$id_vars()
  resData2 <- resData |> dplyr::select(c(id_vars, table$value_vars()) )
  resData2 <- resData2 |> tidyr::gather(key = "variable", value = "value", -dplyr::all_of(id_vars)  )
  resData2 <- resData2 |>  tidyr::unite("temp", table$isotopeLabel, .data$variable )
  HLData <- resData2 |> tidyr::spread(.data$temp,.data$value)
  invisible(HLData)
}


# Computing protein Intensity summaries ---


