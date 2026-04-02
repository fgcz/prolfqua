# Moderation and p-value adjustment ----

#' Moderate p-values - limma approach
#' @export
#' @family modelling
#' @keywords internal
#'
moderated_p_limma <- function(contrast_df, df = "df", estimate = "diff", robust = FALSE, confint = 0.95) {
  squeezed_var <- prolfqua::squeezeVarRob(contrast_df$sigma^2, df = contrast_df[[df]], robust = robust)

  # pior degrees of freedom are Inf
  if (all(is.infinite(squeezed_var$df.prior))) {
    squeezed_var$df.prior <- mean(contrast_df[[df]]) * nrow(contrast_df) / 10
  }

  squeezed_var <- tibble::as_tibble(squeezed_var)
  squeezed_var <- setNames(squeezed_var, paste0("moderated.", names(squeezed_var)))
  contrast_df <- dplyr::bind_cols(contrast_df, squeezed_var)
  contrast_df <- contrast_df |>
    dplyr::mutate(
      moderated.statistic = .data$statistic * .data$sigma / sqrt(.data$moderated.var.post),
      moderated.df.total = !!sym(df) + .data$moderated.df.prior,
      moderated.p.value = 2 * pt(abs(.data$moderated.statistic), df = .data$moderated.df.total, lower.tail = FALSE)
    )

  conf_quantile <- -qt((1 - confint) / 2, df = contrast_df$moderated.df.total)
  contrast_df$moderated.conf.low <- contrast_df[[estimate]] - conf_quantile * sqrt(contrast_df$moderated.var.post)
  contrast_df$moderated.conf.high <- contrast_df[[estimate]] + conf_quantile * sqrt(contrast_df$moderated.var.post)
  contrast_df <- dplyr::ungroup(contrast_df)

  return(contrast_df)
}

#' Moderate p-value for long table
#' @param contrast_df result of \code{\link{contrasts_linfct}}
#' @param group_by_col colnames with contrast description - default 'lhs'
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#'
#' mod <- sim_build_models_lm()
#' m <- get_complete_model_fit(mod$modelDF)
#' factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
#' factor_levelContrasts <- contrasts_linfct(
#'   mod$modelDF,
#'   factor_contrasts,
#'   subject_Id = "protein_Id",
#'   contrastfun = compute_contrast)
#'
#' mmm <- moderated_p_limma_long(factor_levelContrasts, group_by_col = "lhs")
#'
moderated_p_limma_long <- function(contrast_df, group_by_col = "lhs", estimate = "estimate", robust = FALSE) {
  split_groups <- contrast_df |>
    dplyr::group_by(across(all_of(group_by_col))) |>
    dplyr::group_split()
  moderated_results <- purrr::map_df(split_groups, moderated_p_limma, estimate = estimate, robust = robust)
  return(moderated_results)
}


#' adjust columns
#'
#' @export
#' @param contrast_df data.frame with p-values to adjust
#' @param column name of column containing p-values
#' @param group_by_col column(s) to group by before adjusting (e.g. contrast), or NULL for no grouping
#' @param newname name of the new column with adjusted p-values
#' @examples
#'
#' bb <- c(runif(1000), rexp(1500,rate=5))
#' length(bb)
#' bb <- bb[bb < 1]
#' length(bb)
#' bb <- bb[1:2000]
#' hist(bb)
#' data <- data.frame(contrast = rep(LETTERS[1:5],400), p.value = bb)
#'
#' dataX <- adjust_p_values(data)
#' Adata <- dataX |> dplyr::filter(contrast == "A")
#' stopifnot(all.equal(Adata$FDR, p.adjust(Adata$p.value, method="BH")))
#' data2 <- adjust_p_values(data, group_by_col = NULL)
#' stopifnot(all.equal(data2$FDR, p.adjust(data2$p.value, method="BH")))
#'
#'
adjust_p_values <- function(
  contrast_df,
  column = "p.value",
  group_by_col = "contrast",
  newname = "FDR"
) {
  grouped_df <- contrast_df |>
    dplyr::group_by(across(all_of(group_by_col)))
  adjusted_df <- dplyr::mutate(grouped_df, !!newname := p.adjust(!!sym(column), method = "BH"))
  return(adjusted_df)
}


# ROPECA ----

#' p-value of protein from p.value of the median fold change peptide.
#' @param max.n limit number of peptides per protein.
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' plot(get_p_values_pbeta(0.1,1:10,10), ylim=c(0,0.1))
#' plot(get_p_values_pbeta(0.1,1:10,3), ylim=c(0,0.1))
#' plot(get_p_values_pbeta(0.3,1:30, 3), ylim=c(0,0.1))
#' abline(h=.05,col = 2)
#' plot(seq(0.0,1.0,length=30),get_p_values_pbeta(seq(0.0,1.0,length=30),rep(10,30)))
#' abline(0,1)
#' plot(seq(0.0,1.0,length=30),get_p_values_pbeta(seq(0.0,1.0,length=30),rep(10,30),3))
#' abline(0,1)
#' testthat::expect_equal(get_p_values_pbeta(0.3,10, 3),0.216, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(0,10, 3),0, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(1,10, 3),1, tolerance = 1e-4)
#' testthat::expect_equal(get_p_values_pbeta(1,10, 3),get_p_values_pbeta(1,3, 10), tolerance = 1e-4)
#'
get_p_values_pbeta <- function(median.p.value, n.obs, max.n = 10) {
  n.obs <- pmin(n.obs, max.n)

  shape1 <- (n.obs / 2 + 0.5)
  shape2 <- (n.obs - (n.obs / 2 + 0.5) + 1)

  stopifnot(shape1 == shape2)
  res.p.value <- pbeta(median.p.value, shape1 = shape1, shape2 = shape2)
  return(res.p.value)
}


#' compute protein level fold changes and p.values (using beta distribution)
#' takes p-value of the scaled p-value
#'
#' @param contrasts_data data frame
#' @param contrast name of column with contrast identifier
#' @param subject_Id name of column with typically protein Id
#' @param estimate name of column with effect size estimate
#' @param statistic statistic name of column with statistic (typically t-statistics)
#' @param p.value name of column with moderated.p.value
#' @param max.n used to limit the number of peptides in probablity computation.
#' @export
#' @family modelling
#' @keywords internal
#' @return data.frame with columns
#'
#'
#' @examples
#'
#' set.seed(10)
#' nrPep <- 10000
#' nrProtein <- 800
#' p.value <- runif(nrPep)
#' estimate <- rnorm(nrPep)
#' avgAbd <- runif(nrPep)
#' protein_Id <- sample(1:800, size = nrPep,
#'   replace = TRUE, prob = dexp(seq(0,5,length = 800)))
#'
#' plot(table(table(protein_Id)))
#'
#' testdata <- data.frame(contrast = "contrast1",
#'   protein_Id = protein_Id,
#'   estimate = estimate,
#'   pseudo_estimate = estimate,
#'   p.value = p.value,
#'   avgAbd = avgAbd )
#'
#' xx30 <- summary_ROPECA_median_p.scaled(testdata,
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value = "p.value",
#'                                     max.n = 30)
#'
#' xx2 <- summary_ROPECA_median_p.scaled(testdata,
#'                                     subject_Id = "protein_Id",
#'                                     estimate = "estimate",
#'                                     p.value = "p.value",
#'                                     max.n = 1)
#'
#' testthat::expect_equal(mad(xx2$estimate, na.rm = TRUE),0.384409, tolerance = 1e-4)
#' testthat::expect_equal(median(xx2$estimate), -0.006874857, tolerance = 1e-4)
#' testthat::expect_equal(xx2$beta.based.significance[1],0.819, tolerance = 1e-3)
#' testthat::expect_equal(xx2$beta.based.significance[2],0.9234362,tolerance = 1e-3)
#'
#' # Uniform distribution
#' hist(testdata$p.value)
#' hist(xx30$median.p.scaled, breaks = 20)
#' hist(xx2$median.p.scaled, breaks = 20)
#' # shows that beta.based.significance has NO uniform distribution
#' # although H0 is true for all cases.
#'
#' hist(xx30$beta.based.significance, breaks = 20)
#' hist(xx2$beta.based.significance, breaks = 20)
#'
#' hist(xx2$median.p.value, breaks = 20)
#' hist(xx2$beta.based.significance, breaks = 20)
#' hist(estimate)
#'
summary_ROPECA_median_p.scaled <- function(
  contrasts_data,
  contrast = "contrast",
  subject_Id = "protein_Id",
  estimate = "diff",
  statistic = "statistic",
  p.value = "moderated.p.value",
  max.n = 10
) {
  nrpeps_per_prot <- contrasts_data |>
    group_by(across(all_of(c(subject_Id, contrast)))) |>
    dplyr::summarize(n = dplyr::n())

  contrasts_data <- contrasts_data |>
    dplyr::mutate(
      scaled.p = ifelse(!!sym(estimate) > 0, 1 - !!sym(p.value), !!sym(p.value) - 1)
    )

  summarized.protein <- contrasts_data |>
    group_by(across(all_of(c(subject_Id, contrast)))) |>
    dplyr::summarize(
      n_not_na = n(),
      mad.estimate = mad(!!sym(estimate), na.rm = TRUE),
      estimate = median(!!sym(estimate), na.rm = TRUE),
      statistic = median(!!sym(statistic), na.rm = TRUE),
      median.p.scaled = median(.data$scaled.p, na.rm = TRUE),
      avgAbd = median(.data$avgAbd, na.rm = TRUE)
    )

  summarized.protein <- summarized.protein |>
    dplyr::mutate(median.p.value = 1 - abs(.data$median.p.scaled))

  summarized.protein <- summarized.protein |>
    dplyr::mutate(beta.based.significance = get_p_values_pbeta(.data$median.p.value, .data$n_not_na, max.n = max.n))
  summarized.protein <- summarized.protein |>
    dplyr::mutate(n.beta = pmin(.data$n_not_na, max.n))

  summarized.protein <- dplyr::inner_join(nrpeps_per_prot, summarized.protein, by = c(subject_Id, contrast))

  summarized.protein$isSingular <- FALSE
  # scale it back here.
  return(ungroup(summarized.protein))
}


#' Fishers exact test on a datframe
#' @export
#' @family modelling
#' @keywords internal
#' @examples
#' Nprot <- 1000
#' condA <- 8
#' condB <- 8
#' observedA <- sample(0:8, Nprot, replace = TRUE)
#' observedB <- sample(0:8, Nprot, replace = TRUE)
#' fisher_input <- data.frame(observedA = observedA, observedB = observedB)
#'
#' fisher_input$samplesA <- condA
#' fisher_input$samplesB <- condB
#' proteinID <- unique(stringi::stri_rand_strings(Nprot + 20,5))[1:Nprot]
#' fisher_input$proteinID <- proteinID
#' res <- contrasts_fisher_exact(fisher_input)
#'
contrasts_fisher_exact <- function(
  fisher_input,
  observedA = "observedA",
  observedB = "observedB",
  samplesA = "samplesA",
  samplesB = "samplesB"
) {
  relativeRisk <- function(observedA, observedB, samplesA, samplesB) {
    rr <- (observedA / (observedA + observedB)) / (samplesA / (samplesA + samplesB))
    return(rr)
  }
  odsRatio <- function(observedA, observedB, samplesA, samplesB) {
    rr <- (observedA / observedB) / (samplesA / samplesB)
    return(rr)
  }
  apply_fischer <- function(proteinID, observedA, observedB, samplesA, samplesB) {
    mat <- matrix(c(observedA, samplesA - observedA, observedB, samplesB - observedB), nrow = 2)
    fisher_result <- fisher.test(mat)
    return(data.frame(
      proteinID = proteinID,
      p_value = fisher_result$p.value,
      OdsRatio = (fisher_result$estimate),
      conf.lower = (fisher_result$conf.int[1]),
      conf.higher = (fisher_result$conf.int[2])
    ))
  }

  fisher_input$OdsRatioM <- odsRatio(
    observedA = fisher_input[["observedA"]],
    observedB = fisher_input[["observedB"]],
    samplesA = fisher_input[["samplesA"]],
    samplesB = fisher_input[["samplesB"]]
  )
  fisher_input$relativeRiskM <- relativeRisk(
    observedA = fisher_input[["observedA"]],
    observedB = fisher_input[["observedB"]],
    samplesA = fisher_input[["samplesA"]],
    samplesB = fisher_input[["samplesB"]]
  )

  res <- vector(mode = "list", length(nrow(fisher_input)))

  for (i in seq_len(nrow(fisher_input))) {
    res[[i]] <- apply_fischer(
      fisher_input[["proteinID"]][i],
      fisher_input[["observedA"]][i],
      fisher_input[["observedB"]][i],
      fisher_input[["samplesA"]][i],
      fisher_input[["samplesB"]][i]
    )
  }

  result <- dplyr::bind_rows(res)
  enriched_result <- dplyr::inner_join(fisher_input, result, by = c("proteinID" = "proteinID"))
  return(enriched_result)
}
