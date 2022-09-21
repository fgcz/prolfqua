



# AnalysisConfiguration ----
#' Analysis Configuration
#' @keywords internal
#' @family configuration
#' @export
#'
AnalysisConfiguration <- R6::R6Class("AnalysisConfiguration",
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
#'
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
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
#'
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' analysis <- bb$data
#'
#' config$table$factorKeys()
#' config$table$factorDepth <- 1
#' make_interaction_column(analysis,
#'    config$table$fkeysDepth())
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
R6extractValues <- function(r6class){
  tmp <- sapply(r6class, class)
  slots <- tmp[ !tmp %in% c("environment", "function")]
  res <- list()
  for (i in names(slots)) {
    if ("R6" %in% class(r6class[[i]])) {
      res[[i]]  <- R6extractValues(r6class[[i]])
    }else{
      res[[i]] <- r6class[[i]]
    }
  }
  return(res)
}



#' Extracts columns relevant for a configuration from a data frame
#' and create new columns i sampleName column etc.
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
setup_analysis <- function(data, configuration, cc = TRUE ){
  configuration <- configuration$clone(deep = TRUE)
  table <- configuration$table
  for (i in seq_along(table$hierarchy))
  {
    data <- tidyr::unite(data, UQ(sym(table$hierarchyKeys()[i])), table$hierarchy[[i]],remove = FALSE, sep = configuration$sep)
  }
  data <- dplyr::select(data , -dplyr::one_of(dplyr::setdiff(unlist(table$hierarchy), table$hierarchyKeys() )))

  for (i in seq_along(table$factors))
  {
    if ( length(table$factors[[i]]) > 1) {
      data <- tidyr::unite(data, UQ(sym(table$factorKeys()[i])), table$factors[[i]],remove = FALSE, sep = configuration$sep)
    } else {
      newname <- table$factorKeys()[i]
      data <- dplyr::mutate(data, !!newname := !!sym(table$factors[[i]]))
    }
  }

  sampleName <- table$sampleName

  if (!sampleName  %in% names(data)) {
    message("creating sampleName")

    data <- data |>  tidyr::unite( UQ(sym( sampleName)) , unique(unlist(table$factors)), remove = TRUE , sep = configuration$sep) |>
      dplyr::select(sampleName, table$fileName) |> dplyr::distinct() |>
      dplyr::mutate_at(sampleName, function(x){ x <- make.unique( x, sep = configuration$sep )}) |>
      dplyr::inner_join(data, by = table$fileName)
  } else {
    message(sampleName, " already exists")
  }

  data <- data |> dplyr::select(-dplyr::one_of(dplyr::setdiff(unlist(table$factors), table$factorKeys())))

  # Make implicit NA's explicit

  if (!configuration$table$isotopeLabel %in% colnames(data)) {
    warning("no isotopeLabel column specified in the data, adding column automatically and setting to 'light'.")
    data[[configuration$table$isotopeLabel]] <- "light"
  }

  data <- data |> dplyr::select(c(configuration$table$idVars(),configuration$table$valueVars()))

  if (cc) {
    data <- complete_cases( data , configuration)
  }
  return( data )
}

#' separates hierarchy columns into starting columns
#'
#' @export
#' @param data data.frame
#' @param config AnlalysisConfiguration
#' @family configuration
#' @keywords internal
separate_hierarchy <- function(data, config){
  for (hkey in config$table$hkeysDepth()) {
    if (length(config$table$hierarchy[[hkey]]) == 1 & hkey == config$table$hierarchy[[hkey]]) {
      data <- data
    }else {
      data <- data |> tidyr::separate( hkey, config$table$hierarchy[[hkey]], sep = config$sep, remove = FALSE)
    }
  }
  return(data)
}

#' separates factor columns into starting columns
#' @param data data.frame
#' @param config AnlalysisConfiguration
#' @export
#'
#' @keywords internal
#' @family configuration
separate_factors <- function(data, config) {
  for (fkey in config$table$factorKeys()) {
    data <- data |> tidyr::separate( fkey, config$table$factors[[fkey]], sep = config$sep, remove = FALSE)
  }
  return(data)
}


#' Complete cases
#' @export
#' @param pdata data.frame
#' @param config AnlalysisConfiguration
#' @keywords internal
#' @family configuration
#' @examples
#'
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' xx <- complete_cases(data, config)
#' stopifnot(nrow(data) <= nrow(xx))
#'
complete_cases <- function(pdata, config) {
  message("completing cases")
  fkeys <- c(config$table$fileName,config$table$sampleName, config$table$factorKeys())
  hkeys <- c(config$table$isotopeLabel, config$table$hierarchyKeys())
  res <- tidyr::complete(
    pdata,
    tidyr::nesting(!!!syms(fkeys)),
    tidyr::nesting(!!!syms(hkeys))
  )
  return(res)
}


#' sample subset of data
#' @export
#' @keywords internal
#' @family configuration
#'
sample_subset <- function(size, pdata, config){
  hk <- config$table$hkeysDepth()
  message("Sampling ", size, paste(hk, collapse = "," ) )
  hkdf <- pdata |> select(all_of(hk)) |> distinct() |> sample_n(size = size)
  sdata <- inner_join(hkdf, pdata)
  return(sdata)
}

# Functions - summary ----

# Functions - summarize factors ----

#' table of distinct factors (sample annotation)
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#'
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' xx <- table_factors(data,config )
#' xx
#' xx |> dplyr::group_by(!!!rlang::syms(config$table$factorKeys())) |>
#'  dplyr::summarize(n = dplyr::n())
#'
table_factors <- function(pdata, configuration){
  factorsTab <- pdata |> dplyr::select(c(configuration$table$fileName, configuration$table$sampleName, configuration$table$factorKeys())) |>
    dplyr::distinct() |>
    arrange(!!sym(configuration$table$sampleName))
  return(factorsTab)
}


# Functions - summarize hierarchies

#' Count distinct elements for each level of hierarchy
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#'
#' @export
#' @keywords internal
#' @family summary
#' @examples
#'
#'
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
#' stopifnot(nrow(bb$data) == 25780)
#' config <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' x <- hierarchy_counts(data, config)
#' x$protein_Id
#' data <- data |> dplyr::filter( protein_Id == "XYZ")
#' tmp <- hierarchy_counts(data, config)
hierarchy_counts <- function(pdata, config){
  hierarchy <- config$table$hierarchyKeys()
  res <- pdata |>
    dplyr::group_by_at(config$table$isotopeLabel) |>
    dplyr::summarise_at( hierarchy, n_distinct )

  return(res)
}

#' Count distinct elements for each level of hierarchy per sample
#' @export
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#'
#' @keywords internal
#' @family summary
#' @examples
#'
#'
#'
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
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
  hierarchy <- configuration$table$hierarchyKeys()
  summary <- pdata |> dplyr::filter(!is.na(!!sym(configuration$table$getWorkIntensity() ))) |>
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
        nudgeval <- mean(long$nr) * 0.05
        # TODO(WEW) check potential problem with sampleName
        ggplot2::ggplot(long, aes(x = !!sym(configuration$table$sampleName), y = .data$nr)) +
          ggplot2::geom_bar(stat = "identity", position = "dodge", colour = "black", fill = "white") +
          ggplot2::facet_wrap( ~ key, scales = "free_y", ncol = 1) +
          ggplot2::geom_text(aes(label = .data$nr), nudge_y = nudgeval, angle = 45) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      }
    }
  }
  return(res)
}



#' Summarize hierarchy counts
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
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' summarize_hierarchy(data, configur)
#' summarize_hierarchy(data, configur, factors = character())
#'
#' summarize_hierarchy(data, configur,
#'  hierarchy = configur$table$hkeysDepth() )
#' summarize_hierarchy(data, configur,
#'  hierarchy = NULL, factors = configur$table$fkeysDepth() )
#' configur$table$hierarchyDepth = 1
#' summarize_hierarchy(data, configur,
#'  factors = configur$table$fkeysDepth())
#' configur$table$hierarchyDepth = 2
#' summarize_hierarchy(data, configur)
#' configur$table$hierarchyDepth = 3
#' summarize_hierarchy(data, configur )
#' configur$table$hierarchyDepth = 4
#' summarize_hierarchy(data, configur )
#'
summarize_hierarchy <- function(pdata,
                                config,
                                hierarchy = config$table$hkeysDepth(),
                                factors = character())
{
  all_hierarchy <- c(config$table$isotopeLabel, config$table$hierarchyKeys() )

  precursor <- pdata |> dplyr::select(factors, all_hierarchy) |> dplyr::distinct()
  x3 <- precursor |> dplyr::group_by_at(c(factors, hierarchy)) |>
    dplyr::summarize_at( setdiff(all_hierarchy, hierarchy),
                         list( n = dplyr::n_distinct))
  return(x3)
}


# Functions - Handling isotopes ----

#' spreads isotope label heavy light into two columns
#' @export
#' @keywords internal
#' @family configuration
#' @examples
#'
#' bb <- old2new(prolfqua_data('data_ionstar')$filtered())
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config$clone(deep=TRUE)
#' data <- bb$data
#'
#' x<-spreadValueVarsIsotopeLabel(data,configur)
#' head(x)
#'
#' bb <- prolfqua_data('data_skylineSRM_HL_A')
#' configur <- bb$config_f()
#' data <- bb$analysis(bb$data, configur)
#'
#' bb <- prolfqua_data('data_skylineSRM_HL_A')
#' conf <- bb$config_f()
#' analysis <- bb$analysis(bb$data, bb$config_f())
#' x <- spreadValueVarsIsotopeLabel(analysis, conf)
#' head(x[,5:ncol(x)])
#'
spreadValueVarsIsotopeLabel <- function(resData, config){
  table <- config$table
  idVars <- table$idVars()
  resData2 <- resData |> dplyr::select(c(idVars, table$valueVars()) )
  resData2 <- resData2 |> tidyr::gather(key = "variable", value = "value", -dplyr::all_of(idVars)  )
  resData2 <- resData2 |>  tidyr::unite("temp", table$isotopeLabel, .data$variable )
  HLData <- resData2 |> tidyr::spread(.data$temp,.data$value)
  invisible(HLData)
}


# Computing protein Intensity summaries ---


