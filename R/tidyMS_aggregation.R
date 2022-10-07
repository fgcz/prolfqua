# Functions - Plotting ----
# Plot peptide and fragments
plot_hierarchies_line_default <- function(data,
                                          proteinName,
                                          sample,
                                          intensity,
                                          peptide,
                                          fragment,
                                          factor,
                                          isotopeLabel,
                                          separate = FALSE,
                                          log_y = FALSE,
                                          show.legend = FALSE
) {
  if (length(isotopeLabel)) {
    if (separate) {
      formula <- paste(paste( isotopeLabel, collapse = "+"), "~", paste(factor , collapse = "+"))
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group = fragment,
                                   color = peptide
      ))
    }else{
      formula <- sprintf("~%s",paste(factor, collapse = " + "))
      data <- tidyr::unite(data, "fragment_label", fragment, isotopeLabel, remove = FALSE)
      p <- ggplot(data, aes_string(x = sample,
                                   y = intensity,
                                   group = "fragment_label",
                                   color = peptide
      ))
    }
    p <- p +  geom_point(aes_string(shape = isotopeLabel), show.legend = show.legend) +
      geom_line(aes_string(linetype = isotopeLabel), show.legend = show.legend)
  }else{
    formula <- sprintf("~%s", paste(factor, collapse = " + "))
    p <- ggplot(data, aes_string(x = sample, y = intensity, group = fragment,  color = peptide))
    p <- p +  geom_point(show.legend = show.legend) + geom_line(show.legend = show.legend)
  }

  #p <- ggplot(data, aes_string(x = sample, y = intensity, group = fragment,  color= peptide, linetype = isotopeLabel))
  p <- p + facet_grid(as.formula(formula), scales = "free_x"   )
  p <- p + ggtitle(proteinName)
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")
  if (log_y) {
    p <- p + scale_y_continuous(trans = 'log10')
  }
  return(p)
}

#' Plot peptide intensities of protein as a function of the sample and factor
#'
#' @export
#' @param res data.frame
#' @param proteinName title of plot
#' @param config AnalysisConfiguration
#' @param separate if heavy and light show in one plot or with separate y axis?
#' @family aggregation
#' @family plotting
#'
#' @keywords internal
#' @examples
#'
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config
#' analysis <- bb$data
#'
#' xnested <- analysis |>
#'  dplyr::group_by_at(conf$table$hierarchy_keys_depth()) |> tidyr::nest()
#'
#' prolfqua::plot_hierarchies_line(xnested$data[[1]], xnested$protein_Id[[1]],conf )
#'
#' bb <- prolfqua_data('data_skylineSRM_HL_A')
#' conf <- bb$config_f()
#' analysis <- bb$analysis(bb$data, conf)
#'
#' nest <- analysis |> dplyr::group_by(conf$table$hierarchy_keys_depth()) |> tidyr::nest()
#' prolfqua::plot_hierarchies_line(nest$data[[1]],
#'                                   "DUM",
#'                                   conf,
#'                                   separate = TRUE)
#' prolfqua::plot_hierarchies_line(nest$data[[1]],
#' "DUM",
#' conf,
#' separate = TRUE,
#' show.legend = TRUE)
#'
plot_hierarchies_line <- function(res,
                                  proteinName,
                                  config,
                                  separate = FALSE,
                                  show.legend = FALSE){

  rev_hnames <- config$table$hierarchy_keys(TRUE)
  fragment <- rev_hnames[1]
  peptide <- rev_hnames[1]

  if (length(rev_hnames) > 2) {
    peptide <- rev_hnames[2]
  }
  res <- plot_hierarchies_line_default(
    res,
    proteinName = proteinName,
    sample = config$table$sampleName,
    intensity = config$table$get_response(),
    peptide = peptide,
    fragment = fragment,
    factor = config$table$factor_keys_depth(),
    isotopeLabel = config$table$isotopeLabel,
    separate = separate,
    log_y = !config$table$is_response_transformed,
    show.legend = show.legend
  )
  return(res)
}




#' Generates peptide level plots for all proteins
#'
#' @seealso \code{\link{plot_hierarchies_line}}
#' @export
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @family aggregation
#' @family plotting
#' @keywords internal
#' @examples
#'
#'
#' istar <- prolfqua_data('data_ionstar')$normalized()
#'
#' istar_data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 20))
#' config <-  old2new(istar$config)
#'
#' config$table$is_response_transformed <- FALSE
#' #debug(plot_hierarchies_line_df)
#' res <- plot_hierarchies_line_df(istar_data, config)
#' res[[1]]
#'
#' config$table$is_response_transformed <- TRUE
#' res <- plot_hierarchies_line_df(istar_data, config)
#' res[[1]]
#'
#' istar <- prolfqua_data('data_ionstar')$filtered()
#' istar_data <- istar$data |> dplyr::filter(protein_Id %in% sample(protein_Id, 20))
#' config <-  old2new(istar$config)
#' res <- plot_hierarchies_line_df(istar_data, config)
#' config$table$is_response_transformed
#' res[[1]]
#' config$table$is_response_transformed <- TRUE
#' res <- plot_hierarchies_line_df(istar_data, config)
#' config$table$is_response_transformed
#' res[[1]]
#'
#' #TODO make it work for other hiearachy levels.
#' config$table$hierarchyDepth = 2
#' #res <- plot_hierarchies_line_df(resDataStart, config)
#'
plot_hierarchies_line_df <- function(pdata, config, show.legend = FALSE){
  factor_level <- config$table$factorDepth

  hierarchy_ID <- "hierarchy_ID"
  pdata <- pdata |> tidyr::unite(hierarchy_ID , !!!syms(config$table$hierarchy_keys_depth()), remove = FALSE)

  xnested <- pdata |> dplyr::group_by_at(hierarchy_ID) |> tidyr::nest()

  figs <- xnested |>
    dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                              plot_hierarchies_line,
                              config = config, show.legend = show.legend ) )
  return(figs$plot)
}



#' Add protein estimate to plot of peptide intensities
#' @export
#' @family aggregation
#' @family plotting
#' @keywords internal
#' @examples
#' #todo
plot_hierarchies_add_quantline <- function(p, data, aes_y,  configuration){
  table <- configuration$table
  p + geom_line(data = data,
                aes_string(x = table$sampleName , y = aes_y, group = 1),
                size = 1.3,
                color = "black",
                linetype = "solid") +
    geom_point(data = data,
               aes_string(x = table$sampleName , y = aes_y, group = 1), color = "black", shape = 10)
}


.reestablish_condition <- function(data,
                                   medpolishRes,
                                   config
){
  table <- config$table
  xx <- data |>  dplyr::select(c(table$sampleName,
                                  table$factor_keys(),
                                  table$fileName,
                                  table$isotopeLabel)) |> dplyr::distinct()
  res <- dplyr::inner_join(xx,medpolishRes, by = table$sampleName)
  res
}



#' Median polish estimate of e.g. protein from peptide intensities
#'
#' Compute Tukey's median polish estimate of a protein from peptide or precursor intensities
#'
#' @family aggregation
#' @param x a matrix
#' @param name if TRUE returns the name of the summary column
#' @return data.frame with number of rows equal to number of columns of input matrix.
#' @export
#' @keywords internal
#'
#' @examples
#'
#' medpolish_estimate(name = TRUE)
#' gg <- matrix(runif(20),4,5)
#' rownames(gg) <- paste0("A",1:4)
#' colnames(gg) <- make.names(1:5)
#' gg
#' mx <- medpolish_estimate(gg)
#'
medpolish_estimate <- function(x, name = FALSE, sampleName = "sampleName" ){
  if (name) {
    return("medpolish")
  }
  X <- medpolish(x,na.rm = TRUE, trace.iter = FALSE, maxiter = 10);
  res <- tibble(!! sampleName := names(X$col) , medpolish = X$col + X$overall)
  res
}

.extractInt <- function(pdata, expression, feature, sampleName ){
  pdata <- pdata |>
    dplyr::select_at( c( sampleName,
                         feature,
                         expression) ) |>
    tidyr::spread(key = sampleName , value = expression) |> .ExtractMatrix()
  return(pdata)
}

#' Extract response column of a protein into matrix
#'
#' Used to apply the median polish function working on matrices to a tidy table
#'
#' @export
#' @keywords internal
#' @examples
#' library(dplyr)
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' configur <- bb$config
#' data <- bb$data
#'
#' xnested <- data |>
#'  dplyr::group_by_at( configur$table$hierarchy_keys_depth() ) |>
#'  tidyr::nest()
#' x <- xnested$data[[1]]
#' nn  <- x |> dplyr::select( setdiff(configur$table$hierarchy_keys() ,  configur$table$hierarchy_keys_depth()) ) |>
#'  dplyr::distinct() |> nrow()
#'
#' xx <- response_as_matrix(x,configur)
#' stopifnot(dim(xx)==c(nn,20))
#'
#' # change hierarchyDepth ###################
#' conf <- configur$clone(deep=TRUE)
#' conf$table$hierarchyDepth = 1
#'
#' xnested <- data |>
#'  dplyr::group_by_at(conf$table$hierarchy_keys_depth()) |>
#'  tidyr::nest()
#'
#' x <- xnested$data[[1]]
#' nn  <- x |> dplyr::select( setdiff(configur$table$hierarchy_keys(),  configur$table$hierarchy_keys_depth()) ) |>
#'  dplyr::distinct() |> nrow()
#'
#' xx <- response_as_matrix(x,conf)
#' stopifnot(dim(xx)==c(nn,20))
#'
response_as_matrix <- function(pdata, config ){
  table <- config$table
  .extractInt(pdata, table$get_response(),
              setdiff(table$hierarchy_keys(), table$hierarchy_keys_depth()),
              table$sampleName)
}

#' Median polish estimates of e.g. protein abundances for entire data.frame
#'
#' @seealso \code{\link{medpolish_estimate}}
#' @param pdata data.frame
#' @param expression column name with intensities
#' @param feature column name e.g. peptide ids
#' @param samples column name e.g. sampleName
#' @return data.frame
#' @export
#' @keywords internal
#' @family aggregation
#' @family plotting
#' @examples
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config
#' data <- bb$data
#'
#' conf$table$hierarchyDepth = 1
#' xnested <- data |>
#'   dplyr::group_by_at(conf$table$hierarchy_keys_depth()) |> tidyr::nest()
#'
#' feature <- setdiff(conf$table$hierarchy_keys(),  conf$table$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- medpolish_estimate_df(x,
#'  expression = conf$table$get_response(),
#'   feature = feature,
#'    sampleName = conf$table$sampleName)
#' prolfqua:::.reestablish_condition(x,bb, conf)
#'
medpolish_estimate_df <- function(pdata, expression, feature, sampleName  ){
  bb <- .extractInt(pdata,
    expression = expression,
     feature = feature,
      sampleName =  sampleName)
   medpolish_estimate(bb, sampleName = sampleName)
}


#' Median polish estimates of e.g. protein abundances for entire data.frame
#'
#' @seealso \code{\link{medpolish_estimate_df}}
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config
#' data <- bb$data
#' conf$table$hierarchyDepth = 1
#' xnested <- data |>
#'   dplyr::group_by_at(conf$table$hierarchy_keys_depth()) |> tidyr::nest()
#'
#' feature <- setdiff(conf$table$hierarchy_keys(),  conf$table$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- medpolish_estimate_dfconfig(x,conf)
#' prolfqua:::.reestablish_condition(x,bb, conf)
#'
medpolish_estimate_dfconfig <- function(pdata, config, name=FALSE){
  if (name) {
    return("medpolish")
  }

  feature <- setdiff(config$table$hierarchy_keys(),  config$table$hierarchy_keys_depth())
  res <- medpolish_estimate_df(pdata,
                 expression = config$table$get_response(),
                 feature = feature,
                 sampleName = config$table$sampleName)
  return(res)
}


.rlm_estimate <- function(pdata, expression, feature , samples = "samples", maxIt = 20) {
  data <- pdata |> select_at(c(samples, feature, expression)) |> na.omit()
  ##If there is only one 1 peptide for all samples return expression of that peptide
  expname <- paste0("mean.",expression)

  if (length(unique(data[[feature]])) == 1L) {
    data$lmrob <- data[[expression]]
    data$weights <- 1
    data <- rename(data,  !!expname := !!sym(expression))
    data <- data |> select(-!!sym(feature))
    return(data)
  }

  ## model-matrix breaks on factors with 1 level so make vector of ones (will be intercept)
  if (length(unique(data[[samples]])) == 1L) {
    data <- data |> group_by_at(samples) |>
      summarize(lmrob = mean(!!sym(expression)),
                !!expname := mean(!!sym(expression)), .groups = "drop")
    data$weights <- 1
    return(data)
  }

  ## sum contrast on peptide level so sample effect will be mean over all peptides instead of reference level
  formula <- as.formula(paste0("~ -1 + " , samples, " + ", feature ))
  contr.arg <- list('contr.sum')
  names(contr.arg) <- feature
  X = model.matrix(formula , data = data, contrasts.arg = contr.arg)
  ## MASS::rlm breaks on singular values.
  ## check with base lm if singular values are present.
  ## if so, these coefficients will be zero, remove this column from model matrix
  ## rinse and repeat on reduced model-matrix untill no singular values are present
  y <- data[[expression]]
  repeat {
    fit = .lm.fit(X,y)
    id = fit$coefficients != 0
    X = X[ , id, drop = FALSE]
    if (!any(!id)) break
  }
  ## Last step is always rlm if X > has some columns left
  if (ncol(X) > 0) {
    fit = MASS::rlm(X, y, maxit = maxIt)
    data$residuals <- fit$residuals
    usamples <- unique(data[[samples]])
    coefNames <- paste0(samples, usamples)
    lmrob <- tibble(!!samples := usamples, lmrob = fit$coefficients[coefNames])

    sumdata <- data |>
      select(-!!sym(feature)) |>
      group_by_at(samples) |>
      dplyr::summarize(!!expname := mean(!!sym(expression)),
                       weights = 1 / mean(residuals^2), .groups = "drop")
    if (any(is.infinite(sumdata$weights) | is.na(sumdata$weights) | sumdata$weights > 10e6)) {
      sumdata$weights = 1
    }

    res <- inner_join(sumdata, lmrob, by = samples)
    if (sum(is.na(res[[expname]])) < sum(is.na(res$lmrob))) {
      res$lmrob <- res[[expname]]
    }
  } else {
    sumdata <- data |>
      select(-!!sym(feature)) |>
      group_by_at(samples) |>
      dplyr::summarize(!!expname := mean(!!sym(expression)),
                       lmrob = mean(!!sym(expression)),
                       .groups = "drop")
    sumdata$weights <- 1
    res <- sumdata
  }
  pdata <- pdata |> dplyr::select_at(samples) |> distinct()
  res <- left_join(pdata, res, by = samples)
  return(res)
}



#' Estimate e.g. protein abundance from peptides using MASS:rlm
#'
#' @keywords internal
#' @family aggregation
#' @param pdata data
#' @param expression intensities
#' @param feature e.g. peptideIDs.
#' @param sample e.g. sampleName
#' @importFrom MASS rlm
#' @export
#' @examples
#'
#' xx <- data.frame(expression = rnorm(20,0,10), feature = rep(LETTERS[1:5],4), samples= rep(letters[1:4],5))
#'
#' bb <- rlm_estimate(xx , "expression", "feature", "samples", maxIt = 20)
#'
#' xx2 <- data.frame(log2Area = rnorm(20,0,10), peptide_Id = rep(LETTERS[1:5],4), sampleName = rep(letters[1:4],5))
#' rlm_estimate(xx2, "log2Area", "peptide_Id", "sampleName")
#' rlm_estimate(prolfqua_data('data_checksummarizationrobust87'),"log2Area", "peptide_Id", "sampleName")
#' rlm_estimate(prolfqua_data('data_checksummarizerobust69'),"log2Area", "peptide_Id", "sampleName")
#' res <- vector(100,mode = "list")
#' for (i in seq_len(100)) {
#'   xx3 <- xx2
#'   xx3$log2Area[sample(1:20,sample(1:15,1))] <- NA
#'   res[[i]] <- list(data = xx3, summary = rlm_estimate(xx3, "log2Area", "peptide_Id", "sampleName"))
#' }
#' rlm_estimate(xx2[xx2$peptide_Id == 'A',],"log2Area", "peptide_Id", "sampleName")
#' rlm_estimate(xx2[xx2$sampleName == 'a',],"log2Area", "peptide_Id", "sampleName")
#'
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' bb$config <- old2new(bb$config)
#' stopifnot(nrow(bb$data) == 25780)
#' conf <- bb$config
#' data <- bb$data
#' conf$table$hierarchyDepth = 1
#' xnested <- data |>
#'   dplyr::group_by_at(conf$table$hierarchy_keys_depth()) |> tidyr::nest()
#'
#' feature <- setdiff(conf$table$hierarchy_keys(),  conf$table$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- rlm_estimate(x,
#'  expression = conf$table$get_response(),
#'   feature = feature,
#'    samples = conf$table$sampleName)
#'
#'
rlm_estimate <- function(pdata, expression, feature , samples, maxIt = 20) {
  pdata <- unite(pdata, "feature", all_of(feature))

  res <- .rlm_estimate(pdata, expression, feature = "feature", samples, maxIt = maxIt)
  return(res)
}

#' Estimate protein abundance from peptide abundances using MASS::rlm
#'
#'
#'
#' @seealso \code{\link{rlm_estimate}}
#' @export
#' @param pdata data.frame
#' @param config \code{\link{AnalysisConfiguraton}}
#' @family aggregation
#' @keywords internal
#' @examples
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' conf <- old2new(bb$config)
#' data <- bb$data
#' conf$table$hierarchyDepth = 1
#' xnested <- data |>
#'   dplyr::group_by_at(conf$table$hierarchy_keys_depth()) |> tidyr::nest()
#'
#' feature <- setdiff(conf$table$hierarchy_keys(),  conf$table$hierarchy_keys_depth())
#' x <- xnested$data[[1]]
#' bb <- rlm_estimate_dfconfig(x, conf)
#'
#' prolfqua:::.reestablish_condition(x,bb, conf)
#'
rlm_estimate_dfconfig <- function(pdata, config, name= FALSE){
  if (name) {return("lmrob")}

  feature <- setdiff(config$table$hierarchy_keys(),  config$table$hierarchy_keys_depth())
  rlm_estimate(pdata, expression = config$table$get_response(),
                  feature = feature,
                  samples = config$table$sampleName
                  , maxIt = 20)
}

#' Convert old proflqua configurations (prolfqua 0.4) to new Analysis configurations
#' prolfqua 0.5
#' @param list with data = data.frame and config = AnalysisConfiguration
#' @return  list with data and AnalysisConfiguration
#' @export
#' @examples
#' dd <- prolfqua_data('data_ionstar')$filtered()
#' dd$config <- old2new(dd$config)
#'
old2new <- function(config) {
  ata <- AnalysisTableAnnotation$new()
  ata$factors <- config$table$factors
  ata$factorDepth <- config$table$factorDepth
  ata$workIntensity <- config$table$workIntensity
  ata$hierarchyDepth <- config$table$hierarchyDepth
  ata$hierarchy <- config$table$hierarchy
  ata$isotopeLabel <- config$table$isotopeLabel
  ata$is_response_transformed <- config$table$is_intensity_transformed
  ata$ident_qValue <- config$table$ident_qValue
  ata$sampleName <- config$table$sampleName
  ata$isotopeLabel <- config$table$isotopeLabel
  config <- AnalysisConfiguration$new(ata)
  return(config)
}

#' Aggregates e.g. protein abundances from peptide abundances
#'
#' @seealso \code{\link{medpolish_estimate_dfconfig}} \code{\link{rlm_estimate_dfconfig}}
#' @param func - a function working on a matrix of intensities for each protein.
#' @return returns list with data (data.frame) and config (AnalysisConfiguration)
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' dd <- prolfqua_data('data_ionstar')$filtered()
#' dd$config <- old2new(dd$config)
#'
#' config <- dd$config
#' data <- dd$data
#'
#' data <- prolfqua::transform_work_intensity(data, config, log2)
#' colnames(data)
#' bbMed <- estimate_intensity(data, config, .func = medpolish_estimate_dfconfig)
#'
#' bbRob <- estimate_intensity(data, config, .func = rlm_estimate_dfconfig)
#' names(bbMed$data)
#' names(bbRob$data)
#' length(bbMed$data$medpolish)
#' length(bbRob$data$lmrob)
#' plot(bbMed$data$medpolish, bbRob$data$lmrob, log="xy", pch=".")
#' abline(0,1, col=2)
#' plot(bbMed$data$medpolish[1:100], bbRob$data$lmrob[1:100])
#' abline(0,1)
#'
estimate_intensity <- function(data, config, .func)
{
  makeName <- .func(name = TRUE)
  config <- config$clone(deep = TRUE)

  xnested <- data |> group_by_at(config$table$hierarchy_keys_depth()) |> nest()

  pb <- progress::progress_bar$new(total = nrow(xnested))
  message("starting aggregation")

  res <- vector( mode = "list" , length = nrow(xnested) )
  for (i in seq_len(nrow(xnested))) {
    pb$tick()
    aggr <- .func(xnested$data[[i]], config)
    res[[i]] <- .reestablish_condition(xnested$data[[i]], aggr , config)
  }
  xnested[[makeName]] <- res
  newconfig <- make_reduced_hierarchy_config(config,
                                             workIntensity = .func(name = TRUE),
                                             hierarchy = config$table$hierarchy_keys_depth(names = FALSE))

  unnested <- xnested |>
    dplyr::select_at(c(config$table$hierarchy_keys_depth(), makeName)) |>
    tidyr::unnest(cols = makeName) |>
    dplyr::ungroup()

  return(list(data = unnested, config = newconfig))
}

#' Plot feature data and result of aggregation
#'
#' @param data data.frame before aggregation
#' @param config AnalyisConfiguration
#' @param data_aggr data.frame after aggregation
#' @param config_aggr AnalysisConfiguration of aggregated data
#' @family plotting
#' @family aggregation
#' @keywords internal
#' @export
#' @examples
#'
#' dd <- prolfqua_data('data_ionstar')$filtered()
#'
#' config <- old2new(dd$config)
#' data <- dd$data
#'
#' data <- prolfqua::transform_work_intensity(data, config, log2)
#' bbMed <- estimate_intensity(data, config, .func = medpolish_estimate_dfconfig)
#' tmpMed <- plot_estimate(data, config, bbMed$data, bbMed$config)
#' stopifnot("ggplot" %in% class(tmpMed$plots[[1]]))
#' stopifnot("ggplot" %in% class(tmpMed$plots[[2]]))
#' tmpMed$plots[[3]]
#'
#' bbRob <- estimate_intensity(data, config, .func = rlm_estimate_dfconfig)
#' tmpRob <- plot_estimate(data, config, bbRob$data, bbRob$config)
#' stopifnot("ggplot" %in% class(tmpRob$plots[[1]]))
#' stopifnot("ggplot" %in% class(tmpRob$plots[[2]]))
#' tmpRob$plots[[3]]
#'
#'
plot_estimate <- function(data, config, data_aggr, config_reduced, show.legend= FALSE ){
  hierarchy_ID <- "hierarchy_ID"
  xnested <- data |> group_by(!!!syms(config$table$hierarchy_keys_depth())) |> nest()
  xnested <- xnested |> tidyr::unite(hierarchy_ID , !!!syms(config$table$hierarchy_keys_depth()))
  xnested_aggr <- data_aggr |> group_by(!!!syms(config_reduced$table$hierarchy_keys_depth())) |> nest_by(.key = "other")
  xnested_aggr <- xnested_aggr |> tidyr::unite(hierarchy_ID , !!!syms(config$table$hierarchy_keys_depth()))
  xnested_all <- inner_join(xnested, xnested_aggr , by = hierarchy_ID )


  plots <- vector(mode = "list", length = nrow(xnested_all))

  pb <- progress::progress_bar$new(total = nrow(xnested_all))
  for (i in seq_len(nrow(xnested_all))) {
    p1 <- plot_hierarchies_line(xnested_all$data[[i]],
                                xnested_all[[hierarchy_ID]][i],
                                config = config, show.legend = show.legend)
    p2 <- plot_hierarchies_add_quantline(p1,
                                         xnested_all$other[[i]],
                                         config_reduced$table$get_response(),
                                         config)
    plots[[i]] <- p2
    pb$tick()
  }
  xnested_all$plots <- plots

  return(xnested_all)
}


#' Aggregates top N intensities
#'
#' run \link{rank_peptide_by_intensity} first
#' @param pdata data.frame
#' @param config AnalysisConfiguration
#' @param func function to use for aggregation
#' @param N default 3 top intensities.
#' @return list with data and new reduced configuration (config)
#' @family aggregation
#' @export
#' @keywords internal
#' @examples
#'
#' dd <- prolfqua_data('data_ionstar')$filtered()
#' config <- old2new(dd$config)
#' res <- dd$data
#' ranked <- rank_peptide_by_intensity(res,config)
#'
#' mean_f <- function(x, name = FALSE){
#'  if(name){return("mean")};mean(x, na.rm=TRUE)
#'  }
#' sum_f <- function(x, name =FALSE){
#'  if(name){return("sum")};sum(x, na.rm = TRUE)
#'  }
#'
#' resTOPN <- aggregate_intensity_topN(
#'  ranked,
#'  config,
#'  .func = mean_f,
#'   N=3)
#'
#' print(dim(resTOPN$data))
#' # stopifnot(dim(resTOPN$data) == c(3260, 8))
#' stopifnot( names(resTOPN) %in% c("data", "config") )
#' config$table$get_response()
#' #debug(plot_estimate)
#' tmpRob <- plot_estimate(ranked,
#'  config,
#'  resTOPN$data,
#'  resTOPN$config,
#'  show.legend=TRUE)
#' stopifnot( "ggplot" %in% class(tmpRob$plots[[4]]) )
#'
aggregate_intensity_topN <- function(pdata , config, .func, N = 3){

  xcall <- as.list( match.call() )
  newcol <- make.names(paste0("srm_",.func(name = TRUE),"_",xcall$N))

  topInt <-
    pdata |> dplyr::filter( !!sym("srm_meanIntRank")  <= N )

  topInt <- topInt |>
    dplyr::group_by_at(c( config$table$hierarchy_keys_depth(),
                          config$table$sampleName,
                          config$table$fileName,
                          config$table$isotopeLabel,
                          config$table$factor_keys() ))
  sumTopInt <- topInt |>
    dplyr::summarize( !!newcol := .func(!!sym(config$table$get_response())),
                      ident_qValue = min(!!sym(config$table$ident_qValue)), .groups = "drop")

  newconfig <- make_reduced_hierarchy_config(
    config,
    workIntensity = newcol,
    hierarchy = config$table$hierarchy[seq_len(config$table$hierarchyDepth)])
  return(list(data = sumTopInt, config = newconfig))
}



#' Summarizes the intensities within hierarchy
#'
#' @param func - a function working on a matrix of intensities for each protein.
#' @return retuns function object
#' @keywords internal
#'
#' @family aggregation
#' @family deprecated
#' @export
#' @examples
#'
#'
#' bb <- prolfqua_data('data_ionstar')$filtered()
#' stopifnot(nrow(bb$data) == 25780)
#' config <- old2new(bb$config$clone(deep = TRUE))
#' data <- bb$data
#' x <- intensity_summary_by_hkeys(data, config, func = medpolish_estimate)
#'
#' res <- x("unnest")
#' res$data |> dim()
#'
#' config <- old2new(bb$config$clone(deep = TRUE))
#' config$table$hierarchyDepth <- 1
#'
#' x <- intensity_summary_by_hkeys(data, config, func = medpolish_estimate)
#' dd <- x(value = "plot")
#' stopifnot(nrow(dd) == length(unique(bb$data$protein_Id)))
#'
#' dd$plot[[2]]
#'
#' # example how to add peptide count information
#' tmp <- summarize_hierarchy(data, config)
#' tmp <- dplyr::inner_join(tmp, x("wide")$data, by = config$table$hierarchy_keys_depth())
#'
intensity_summary_by_hkeys <- function(data, config, func)
{
  x <- as.list( match.call() )
  makeName <- make.names(as.character(x$func))
  config <- config$clone(deep = TRUE)

  xnested <- data |> group_by_at(config$table$hierarchy_keys_depth()) |> nest()

  pb <- progress::progress_bar$new(total = 3 * nrow(xnested))
  message("starting aggregation")

  xnested <- xnested |>
    dplyr::mutate(spreadMatrix = map(data, function(x,config){pb$tick(); response_as_matrix(x, config)}, config))

  #xnested <- xnested |>
  #  dplyr::mutate(!!makeName := map( .data$spreadMatrix , function(x){pb$tick(); func(x)}))
  res <- vector(mode="list", length(nrow(xnested)))
  tmp <- function(x){pb$tick(); func(x)}
  for(i in seq_len(nrow(xnested))){
    res[[i]] <- tmp(xnested$spreadMatrix[[i]])
  }
  xnested[[makeName]] = res
  xnested <- xnested |>
    dplyr::mutate(!!makeName := map2(data, !!sym(makeName), function(x, y, config){pb$tick(); .reestablish_condition(x,y, config) }, config ))


  res_fun <- function(value = c("nested", "unnest", "wide", "plot"), DEBUG = FALSE){

    value <- match.arg(value)
    if (DEBUG) {
      return(list(config = config, value = value, xnested = xnested  ))
    }

    newconfig <- make_reduced_hierarchy_config(config,
                                               workIntensity = func(name = TRUE),
                                               hierarchy = config$table$hierarchy_keys_depth(names = FALSE))

    if (value == "nested") {
      return(list(xnested = xnested, config = newconfig))
    }else if (value == "unnest" || value == "wide") {
      unnested <- xnested |>
        dplyr::select(config$table$hierarchy_keys_depth(), makeName) |>
        tidyr::unnest(cols = c(medpolish_estimate)) |>
        dplyr::ungroup()

      if (value == "wide") {
        wide <- prolfqua::tidy_to_wide_config(unnested, newconfig)
        wide$config <- newconfig
        return(wide)
      }
      return(list(data = unnested, config = newconfig))
    }else if (value == "plot") {
      hierarchy_ID <- "hierarchy_ID"
      xnested <- xnested |> tidyr::unite(hierarchy_ID , !!!syms(config$table$hierarchy_keys_depth()))
      figs <- xnested |>
        dplyr::mutate(plot = map2(data, !!sym(hierarchy_ID) ,
                                  plot_hierarchies_line, config = config ))

      figs <- figs |>
        dplyr::mutate(plot = map2(plot, !!sym(makeName) ,
                                  plot_hierarchies_add_quantline, func(name = TRUE), config ))
      return(figs)
    }
  }
  return(res_fun)
}

#' median polish from normalized peptide intensities
#' @export
#' @keywords internal
#' @family aggregation
#' @family deprecated
#' @examples
#'
#' istar <- prolfqua_data('data_ionstar')$normalized()
#' istar$config <- old2new(istar$config)
#' istar_data <- istar$data
#' res <- medpolish_protein_estimates(istar_data,
#' istar$config )
#'
#' dr <- res("unnest")$data
#' stopifnot(nrow(dr) ==
#' length(unique(istar$data$protein_Id )) * length(unique(istar$data$raw.file)))
#'
medpolish_protein_estimates <- function(data, config){
  protintensity <- prolfqua::intensity_summary_by_hkeys(data ,
                                                          config,
                                                          medpolish_estimate)
  return(protintensity)
}




