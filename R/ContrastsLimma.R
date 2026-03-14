# strategy_limma -----

#' Create limma modelling strategy
#'
#' Analogous to \code{\link{strategy_lm}} but for limma's matrix-based pipeline.
#' Returns a list consumed by \code{\link{build_model_limma}}.
#'
#' @param modelstr model formula as string (e.g. "abundance ~ group_")
#' @param model_name name of model
#' @param trend logical, passed to \code{\link[limma]{eBayes}}
#' @param robust logical, passed to \code{\link[limma]{eBayes}}
#' @export
#' @family modelling
#' @examples
#' strat <- strategy_limma("abundance ~ group_")
#' strat$formula
#' strat$model_name
strategy_limma <- function(modelstr, model_name = "limma", trend = FALSE, robust = FALSE) {
  formula <- as.formula(modelstr)
  list(
    formula = formula,
    model_name = model_name,
    trend = trend,
    robust = robust
  )
}


# build_model_limma -----

#' Build limma model from LFQData
#'
#' Analogous to \code{\link{build_model}} but uses limma's matrix-based pipeline.
#' Takes an LFQData object and a strategy from \code{\link{strategy_limma}},
#' pivots data to wide format, fits with \code{\link[limma]{lmFit}}, and returns
#' a \code{\link{ModelLimma}} object.
#'
#' @param lfqdata an \code{\link{LFQData}} object
#' @param strategy output of \code{\link{strategy_limma}}
#' @param modelName name of model (default from strategy)
#' @return a \code{\link{ModelLimma}} object
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lProt <- LFQData$new(istar$data, istar$config)
#' lProt$rename_response("transformedIntensity")
#'
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod_limma <- build_model_limma(lProt, strat)
#' mod_limma$get_coefficients()
#' mod_limma$get_anova()
#'
build_model_limma <- function(lfqdata, strategy, modelName = strategy$model_name) {
  wide <- lfqdata$to_wide(as.matrix = TRUE)
  expr_matrix <- wide$data # rows = proteins, cols = samples
  annotation <- wide$annotation
  subject_Id <- lfqdata$config$hierarchy_keys()
  rowdata <- wide$rowdata |> dplyr::select(dplyr::all_of(subject_Id))
  if (anyDuplicated(rowdata) && !is.null(lfqdata$config$isotopeLabel)) {
    rowdata <- wide$rowdata |> dplyr::select(dplyr::all_of(unique(c(subject_Id, lfqdata$config$isotopeLabel))))
    subject_Id <- colnames(rowdata)
  }

  # Use only RHS of formula for design matrix (response is the expression matrix)
  rhs_formula <- formula(delete.response(terms(strategy$formula)))
  design <- model.matrix(rhs_formula, data = annotation)

  fit <- limma::lmFit(expr_matrix, design)

  # Fit a dummy lm on one protein's complete data for linfct extraction
  # Pick the first protein with all non-NA values (= complete row)
  complete_rows <- which(rowSums(is.na(expr_matrix)) == 0)
  if (length(complete_rows) == 0) {
    # fallback: pick the row with fewest NAs
    complete_rows <- which.min(rowSums(is.na(expr_matrix)))
  }
  idx <- complete_rows[1]

  dummy_data <- annotation
  dummy_data$.response <- as.numeric(expr_matrix[idx, ])
  dummy_formula <- update(rhs_formula, .response ~ .)
  dummy_model <- lm(dummy_formula, data = dummy_data)

  ModelLimma$new(
    fit = fit,
    design = design,
    formula = strategy$formula,
    subject_Id = subject_Id,
    modelName = modelName,
    rowdata = rowdata,
    trend = strategy$trend,
    robust = strategy$robust,
    dummy_model = dummy_model,
    p.adjust = prolfqua::adjust_p_values
  )
}


# ModelLimma -----

#' R6 class representing a limma modelling result
#'
#' Same API as \code{\link{Model}}: \code{get_anova()}, \code{get_coefficients()},
#' \code{coef_histogram()}, \code{coef_volcano()}, \code{anova_histogram()}.
#'
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lProt <- LFQData$new(istar$data, istar$config)
#' lProt$rename_response("transformedIntensity")
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod <- build_model_limma(lProt, strat)
#'
#' coeffs <- mod$get_coefficients()
#' head(coeffs)
#' anova_tbl <- mod$get_anova()
#' head(anova_tbl)
#' mod$coef_histogram()
#' mod$coef_volcano()
#' mod$anova_histogram()
#'
ModelLimma <- R6::R6Class(
  "ModelLimma",
  inherit = ModelInterface,
  public = list(
    #' @field fit limma MArrayLM object from lmFit
    fit = NULL,
    #' @field design design matrix
    design = NULL,
    #' @field formula model formula
    formula = NULL,
    #' @field subject_Id protein ID column name(s)
    subject_Id = character(),
    #' @field modelName model name
    modelName = character(),
    #' @field rowdata protein ID mapping from to_wide()$rowdata
    rowdata = NULL,
    #' @field trend passed to eBayes
    trend = FALSE,
    #' @field robust passed to eBayes
    robust = FALSE,
    #' @field dummy_model one fitted lm for linfct extraction
    dummy_model = NULL,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize ModelLimma
    #' @param fit limma MArrayLM from lmFit
    #' @param design design matrix
    #' @param formula model formula
    #' @param subject_Id protein ID column name(s)
    #' @param modelName model name
    #' @param rowdata protein ID mapping
    #' @param trend passed to eBayes
    #' @param robust passed to eBayes
    #' @param dummy_model one fitted lm for linfct extraction
    #' @param p.adjust function to adjust p-values
    initialize = function(
      fit,
      design,
      formula,
      subject_Id,
      modelName,
      rowdata,
      trend = FALSE,
      robust = FALSE,
      dummy_model = NULL,
      p.adjust = prolfqua::adjust_p_values
    ) {
      self$fit <- fit
      self$design <- design
      self$formula <- formula
      self$subject_Id <- subject_Id
      self$modelName <- modelName
      self$rowdata <- rowdata
      self$trend <- trend
      self$robust <- robust
      self$dummy_model <- dummy_model
      self$p.adjust <- p.adjust
    },
    #' @description
    #' return model coefficient table in long format
    #' @return data.frame
    get_coefficients = function() {
      fit_eb <- limma::eBayes(self$fit, trend = self$trend, robust = self$robust)

      coef_mat <- fit_eb$coefficients
      se_mat <- fit_eb$stdev.unscaled * fit_eb$s2.post^0.5
      t_mat <- fit_eb$t
      p_mat <- fit_eb$p.value

      ncoef <- ncol(coef_mat)
      coef_names <- colnames(coef_mat)
      ngenes <- nrow(coef_mat)

      res_list <- vector("list", ncoef)
      for (j in seq_len(ncoef)) {
        df_j <- self$rowdata
        df_j$factor <- coef_names[j]
        df_j$Estimate <- coef_mat[, j]
        df_j$Std..Error <- se_mat[, j]
        df_j$t.value <- t_mat[, j]
        df_j$Pr...t.. <- p_mat[, j]
        res_list[[j]] <- df_j
      }
      result <- dplyr::bind_rows(res_list)
      return(result)
    },
    #' @description
    #' return anova table (F-test per protein across all non-intercept coefficients)
    #' @return data.frame
    get_anova = function() {
      fit_eb <- limma::eBayes(self$fit, trend = self$trend, robust = self$robust)

      coef_names <- colnames(fit_eb$coefficients)
      non_intercept <- which(coef_names != "(Intercept)")

      if (length(non_intercept) == 0) {
        warning("No non-intercept coefficients for ANOVA.")
        return(data.frame())
      }

      # Use limma topTable with all non-intercept coefficients for F-test
      tt <- limma::topTable(fit_eb, coef = non_intercept, number = Inf, sort.by = "none")

      result <- self$rowdata
      result$F.value <- tt$F
      result$p.value <- tt$P.Value

      # create a factor column describing the tested term
      terms_tested <- paste(coef_names[non_intercept], collapse = "+")
      result$factor <- terms_tested

      result <- self$p.adjust(result, column = "p.value", group_by_col = "factor")
      return(dplyr::ungroup(result))
    },
    #' @description
    #' histogram of model coefficient p-values
    coef_histogram = function() {
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      fname <- paste0("Coef_Histogram_", self$modelName, ".pdf")
      p <- ggplot(data = Model_Coeff, aes(x = .data$Pr...t.., group = .data$factor)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = p, name = fname))
    },
    #' @description
    #' volcano plot of non-intercept coefficients
    coef_volcano = function() {
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      fname <- paste0("Coef_VolcanoPlot_", self$modelName, ".pdf")
      p <- Model_Coeff |>
        dplyr::filter(.data$factor != "(Intercept)") |>
        prolfqua::multigroup_volcano(
          effect = "Estimate",
          significance = "Pr...t..",
          contrast = "factor",
          label = "subject_Id",
          xintercept = c(-1, 1),
          colour = NULL
        )
      return(list(plot = p, name = fname))
    },
    #' @description
    #' pairs plot of coefficients
    coef_pairs = function() {
      Model_Coeff <- self$get_coefficients()
      Model_Coeff <- tidyr::unite(Model_Coeff, "subject_Id", self$subject_Id)
      forPairs <- Model_Coeff |>
        dplyr::select(all_of(c("subject_Id", "factor", "Estimate"))) |>
        tidyr::pivot_wider(names_from = "factor", values_from = "Estimate")
      fname <- paste0("Coef_Pairsplot_", self$modelName, ".pdf")
      return(list(plot = forPairs, name = fname))
    },
    #' @description
    #' histogram of ANOVA F-test p-values
    #' @param what show either "p.value" or "FDR"
    anova_histogram = function(what = c("p.value", "FDR")) {
      what <- match.arg(what)
      Model_Anova <- self$get_anova()
      fname <- paste0("Anova_p.values_", self$modelName, ".pdf")
      p <- Model_Anova |>
        ggplot(aes(x = !!sym(what), group = .data$factor)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05)) +
        facet_wrap(~factor)
      return(list(plot = p, name = fname))
    }
  )
)


# ContrastsLimma -----

#' Limma-based contrasts (direct limma pipeline)
#'
#' Uses limma's \code{contrasts.fit} + \code{eBayes} pipeline directly,
#' rather than fitting per-protein lm models and then moderating.
#' Inherits from \code{\link{ContrastsInterface}} with the same API as
#' \code{\link{Contrasts}}.
#'
#' @export
#' @family modelling
#' @examples
#'
#' istar <- sim_lfq_data_protein_config(Nprot = 50)
#' lProt <- LFQData$new(istar$data, istar$config)
#' lProt$rename_response("transformedIntensity")
#'
#' strat <- strategy_limma("transformedIntensity ~ group_")
#' mod_limma <- build_model_limma(lProt, strat)
#'
#' Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' contr_limma <- ContrastsLimma$new(mod_limma, Contr)
#' res <- contr_limma$get_contrasts()
#' head(res)
#' stopifnot(all(c("diff", "FDR", "p.value", "statistic") %in% colnames(res)))
#'
#' # Compare with prolfqua's own pipeline
#' modelFunction <- strategy_lm("transformedIntensity ~ group_")
#' mod <- build_model(lProt, modelFunction)
#' contr_prolfqua <- Contrasts$new(mod, Contr)
#' res_prolfqua <- contr_prolfqua$get_contrasts()
#'
#' # fold changes should be very similar
#' merged <- dplyr::inner_join(
#'   dplyr::select(res, protein_Id, diff_limma = diff),
#'   dplyr::select(res_prolfqua, protein_Id, diff_prolfqua = diff),
#'   by = "protein_Id")
#' stopifnot(cor(merged$diff_limma, merged$diff_prolfqua, use = "complete.obs") > 0.99)
#'
#' # Plotter works
#' pl <- contr_limma$get_Plotter()
#'
#' # to_wide works
#' wide <- contr_limma$to_wide()
#' head(wide)
#'
#' # merge_contrasts_results works
#' Contr2 <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' csi <- ContrastsMissing$new(lProt, contrasts = Contr2)
#' merged_res <- merge_contrasts_results(contr_limma, csi)
#'
ContrastsLimma <- R6::R6Class(
  "ContrastsLimma",
  inherit = ContrastsInterface,
  public = list(
    #' @field model ModelLimma object
    model = NULL,
    #' @field contrasts named character vector of contrasts
    contrasts = character(),
    #' @field modelName model name
    modelName = character(),
    #' @field subject_Id columns with subject_Id (proteinID)
    subject_Id = character(),
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @field contrast_result cached contrast results
    contrast_result = NULL,
    #' @description
    #' initialize ContrastsLimma
    #' @param model a \code{\link{ModelLimma}} object
    #' @param contrasts named character vector of contrasts
    #' @param p.adjust function to adjust p-values
    #' @param modelName name of the contrast method
    initialize = function(model, contrasts, p.adjust = prolfqua::adjust_p_values, modelName = "limma") {
      self$model <- model
      self$contrasts <- contrasts
      self$subject_Id <- model$subject_Id
      self$modelName <- modelName
      self$p.adjust <- p.adjust
    },
    #' @description
    #' get both sides of contrasts
    get_contrast_sides = function() {
      tt <- self$contrasts[grep("-", self$contrasts)]
      tt <- tibble(contrast = names(tt), rhs = tt)
      tt <- tt |>
        mutate(rhs = gsub("[` ]", "", rhs)) |>
        tidyr::separate(rhs, c("group_1", "group_2"), sep = "-")
      return(tt)
    },
    #' @description
    #' get linear functions from contrasts
    #' @param global ignored (for API compatibility)
    #' @param avg logical, also compute avgAbd linfct
    get_linfct = function(global = TRUE, avg = TRUE) {
      linfct <- .linfct(self$model$dummy_model, self$contrasts, avg = avg)
      return(linfct)
    },
    #' @description
    #' get table with contrast estimates via limma pipeline
    #' @param all should all columns be returned (default FALSE)
    #' @return data.frame with contrasts
    get_contrasts = function(all = FALSE) {
      if (!is.null(self$contrast_result)) {
        return(self$contrast_result)
      }

      # Build the contrast matrix via linfct_from_model + linfct_matrix_contrasts
      linfct <- linfct_from_model(self$model$dummy_model, as_list = FALSE)
      linfct <- unique(linfct) # needed for single factor models

      # Add avg contrasts for avgAbd computation
      namtmp <- paste0("avg_", names(self$contrasts))
      cntr_avg <- paste0("(", gsub(" - ", " + ", self$contrasts), ")/2")
      names(cntr_avg) <- namtmp
      all_contrasts <- c(self$contrasts, cntr_avg)

      linfct_A <- linfct_matrix_contrasts(linfct, all_contrasts)
      # linfct_A: rows = contrast names, cols = model coefficients

      # Split into difference contrasts and avg contrasts
      diff_names <- names(self$contrasts)
      avg_names <- namtmp

      # Transpose for limma: rows = coefficients, cols = contrasts
      contrast_matrix <- t(linfct_A[diff_names, , drop = FALSE])

      # limma pipeline: contrasts.fit + eBayes
      fit2 <- limma::contrasts.fit(self$model$fit, contrast_matrix)
      fit2 <- limma::eBayes(fit2, trend = self$model$trend, robust = self$model$robust)

      # Extract results per contrast
      res_list <- vector("list", length(diff_names))
      for (i in seq_along(diff_names)) {
        tt <- limma::topTable(fit2, coef = i, number = Inf, sort.by = "none", confint = TRUE)

        df_i <- self$model$rowdata
        df_i$contrast <- diff_names[i]
        df_i$diff <- tt$logFC
        df_i$std.error <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, i]
        df_i$statistic <- tt$t
        df_i$p.value <- tt$P.Value
        df_i$sigma <- sqrt(fit2$s2.post)
        df_i$df <- fit2$df.total
        df_i$conf.low <- tt$CI.L
        df_i$conf.high <- tt$CI.R

        res_list[[i]] <- df_i
      }
      contrast_result <- dplyr::bind_rows(res_list)

      # Compute avgAbd from the avg linfct applied to the fit coefficients
      avg_matrix <- t(linfct_A[avg_names, , drop = FALSE])
      avg_vals <- self$model$fit$coefficients %*% avg_matrix
      # avg_vals: rows = proteins, cols = avg contrasts

      avg_df_list <- vector("list", length(avg_names))
      for (i in seq_along(avg_names)) {
        df_avg <- self$model$rowdata
        df_avg$contrast <- diff_names[i] # map back to original contrast name
        df_avg$avgAbd <- avg_vals[, i]
        avg_df_list[[i]] <- df_avg
      }
      avg_df <- dplyr::bind_rows(avg_df_list)

      contrast_result <- dplyr::left_join(
        contrast_result,
        avg_df,
        by = c(self$subject_Id, "contrast")
      )

      # Adjust p-values per contrast
      contrast_result <- self$p.adjust(contrast_result, column = "p.value", group_by_col = "contrast", newname = "FDR")
      contrast_result <- contrast_result |> dplyr::relocate("FDR", .after = "diff")
      contrast_result <- dplyr::mutate(contrast_result, modelName = self$modelName, .before = 1)
      contrast_result <- dplyr::ungroup(contrast_result)
      self$contrast_result <- contrast_result

      stopifnot(all(super$column_description()$column_name %in% colnames(contrast_result)))
      return(contrast_result)
    },
    #' @description
    #' return \code{\link{ContrastsPlotter}}
    #' @param FCthreshold fold change threshold to show in plots
    #' @param FDRthreshold FDR threshold to show in plots
    #' @return \code{\link{ContrastsPlotter}}
    get_Plotter = function(FCthreshold = 1, FDRthreshold = 0.1) {
      contrast_result <- self$get_contrasts()
      res <- ContrastsPlotter$new(
        contrast_result,
        subject_Id = self$subject_Id,
        fcthresh = FCthreshold,
        volcano = list(
          list(score = "p.value", thresh = FDRthreshold),
          list(score = "FDR", thresh = FDRthreshold)
        ),
        histogram = list(
          list(score = "p.value", xlim = c(0, 1, 0.05)),
          list(score = "FDR", xlim = c(0, 1, 0.05))
        ),
        score = list(list(score = "statistic", thresh = 5)),
        modelName = "modelName",
        diff = "diff",
        contrast = "contrast"
      )
      return(res)
    },
    #' @description
    #' convert to wide format
    #' @param columns value column default p.value
    #' @return data.frame
    to_wide = function(columns = c("p.value", "FDR", "statistic")) {
      contrast_minimal <- self$get_contrasts()
      contrasts_wide <- pivot_model_contrasts_2_Wide(
        contrast_minimal,
        subject_Id = self$subject_Id,
        columns = c("diff", columns),
        contrast = "contrast"
      )
      return(contrasts_wide)
    }
  )
)
