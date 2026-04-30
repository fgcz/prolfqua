# ContrastsPlotter ----
#' plot contrasts
#'
#' @return An R6 class generator.
#' @export
#' @family modelling
#' @family plotting
#' @examples
#'
#' istar <- sim_lfq_data_peptide_config(Nprot = 20)
#' model_name <- "Model"
#' modelFunction <-
#'   strategy_lmer("abundance  ~ group_ + (1 | peptide_Id) + (1|sample)",
#'    model_name = model_name)
#' pepIntensity <- istar$data
#' config <- istar$config
#'
#' mod <- build_model(
#'   pepIntensity,
#'   modelFunction,
#'   model_name = model_name,
#'   subject_id = config$hierarchy_keys_depth())
#'
#'  Contr <- c("group_A_vs_Ctrl" = "group_A - group_Ctrl",
#'  "group_B_vs_Ctrl" = "group_B - group_Ctrl"
#'  )
#' contrast <- prolfqua::Contrasts$new(mod,
#'   Contr)
#' contr <- contrast$get_contrasts()
#' cp <- ContrastsPlotter$new(contr,
#'   contrast$subject_id,
#'   volcano = list(list(score = "FDR", thresh = 0.1)),
#'   histogram = list(list(score = "p.value", xlim = c(0,1,0.05)),
#'                  list(score = "FDR", xlim = c(0,1,0.05))),
#'   score =list(list(score = "statistic",  thresh = 5))
#'   )
#'
#' stopifnot("plotly" %in% class(cp$volcano_plotly()$FDR))
#' stopifnot("ggplot" %in% class(cp$score_plot(legend=FALSE)$statistic))
#' p <- cp$histogram()
#' stopifnot("ggplot" %in% class(p$FDR))
#' stopifnot("ggplot" %in% class(p$p.value))
#' p <- cp$histogram_estimate()
#' stopifnot("ggplot" %in% class(p))
#' res <- cp$volcano()
#' stopifnot("ggplot" %in% class(res$FDR))
#' respltly <- cp$volcano_plotly()
#' stopifnot("plotly" %in% class(respltly$FDR))
#' stopifnot("ggplot" %in%  class(cp$ma_plot()))
#' stopifnot("plotly" %in%   class(cp$ma_plotly(rank=TRUE)))
#' res  <- cp$barplot_threshold()
#' stopifnot("ggplot" %in% class(res$FDR$plot))
#' stopifnot("ggplot" %in%  class(cp$histogram_diff()))
#'
ContrastsPlotter <- R6::R6Class(
  "ContrastsPlotter",
  public = list(
    #' @field contrast_df data frame with contrasts
    contrast_df = NULL,
    #' @field model_name of column with model name
    model_name = character(),
    #' @field subject_id hierarchy key columns
    subject_id = character(),
    #' @field prefix default Contrasts - used to generate file names
    prefix = "Contrasts",
    #' @field diff column with fold change differences
    diff = "diff",
    #' @field contrast column with contrasts names, default "contrast"
    contrast = "contrast",
    #' @field volcano_spec volcano plot specification
    volcano_spec = NULL,
    #' @field score_spec score plot specification
    score_spec = NULL,
    #' @field histogram_spec plot specification
    histogram_spec = NULL,
    #' @field fcthresh fold change threshold
    fcthresh = 1,
    #' @field avg.abundance name of column containing avg abundance values.
    avg.abundance = character(),
    #' @field protein_annot protein annotation
    protein_annot = NULL,
    #' @description
    #' create Crontrast_Plotter
    #' @param contrast_df frame with contrast data
    #' @param subject_id columns containing subject Identifier
    #' @param volcano which score to plot and which ablines to add.
    #' @param histogram which scores to plot and which range (x) should be shown.
    #' @param score score parameters
    #' @param fcthresh default 1 (log2 FC threshold)
    #' @param modelName name of column with model names
    #' @param diff fold change (difference) diff column
    #' @param contrast contrast column
    #' @param avg.abundance name of column with average abundance
    #' @param protein_annot add protein annotation (optional)
    initialize = function(
      contrast_df,
      subject_id,
      volcano = list(list(score = "FDR", thresh = 0.1)),
      histogram = list(list(score = "p.value", xlim = c(0, 1, 0.05)), list(score = "FDR", xlim = c(0, 1, 0.05))),
      score = list(list(score = "statistic", thresh = NULL)),
      fcthresh = 1,
      modelName = "modelName",
      diff = "diff",
      contrast = "contrast",
      avg.abundance = "avgAbd",
      protein_annot = NULL
    ) {
      self$contrast_df <- tidyr::unite(
        contrast_df,
        "subject_id",
        dplyr::all_of(subject_id),
        sep = "~",
        remove = FALSE
      )

      self$model_name <- modelName
      self$subject_id <- subject_id
      self$diff <- diff
      self$volcano_spec <- volcano
      self$score_spec <- score
      self$histogram_spec <- histogram
      self$fcthresh <- fcthresh
      self$contrast <- contrast
      self$avg.abundance <- avg.abundance
      self$protein_annot <- protein_annot
    },
    #' @description
    #' plot histogram of selected scores (e.g. p-value, FDR, t-statistics)
    histogram = function() {
      res <- list()
      if (length(self$histogram_spec) > 0) {
        for (spec in self$histogram_spec) {
          fig <- private$.histogram(score = spec)
          res[[spec$score]] <- fig
        }
        return(res)
      }
    },
    #' @description
    #' plot histogram of effect size - difference between groups
    #' @param binwidth with of bin in histogram
    histogram_estimate = function(binwidth = 0.05) {
      re <- range(self$contrast_df[[self$diff]], na.rm = TRUE)
      re[1] <- floor(re[1])
      re[2] <- ceiling(re[2])

      fig <- ggplot(self$contrast_df, aes(x = !!sym(self$diff))) +
        geom_histogram(breaks = seq(from = re[1], to = re[2], by = binwidth)) +
        geom_vline(
          aes(xintercept = median(!!sym(self$diff), na.rm = TRUE)), # Ignore NA values for mean
          color = "red",
          linetype = "dashed",
          linewidth = 1
        ) +
        geom_vline(xintercept = 0, col = "green", linewidth = 1) +
        facet_wrap(vars(!!sym(self$contrast)))

      return(fig)
    },
    #' @description
    #' plot histogram of differences (diff) fold change
    #' @param binwidth with of bin in histogram
    histogram_diff = function(binwidth = 0.05) {
      self$histogram_estimate(binwidth = binwidth)
    },
    #' @description
    #' volcano plots (fold change vs FDR)
    #' @param colour column name with color information default modelName
    #' @param legend default TRUE
    #' @param scales default fixed \code{\link[ggplot2]{facet_wrap}}, scales argument
    #' @param min_score replace p.values or FDR's smaller then min_score with min_score (default 0.0001).
    volcano = function(colour, legend = TRUE, scales = c("fixed", "free", "free_x", "free_y"), min_score = 0.0001) {
      if (missing(colour)) {
        colour <- self$model_name
      }
      scales <- match.arg(scales)
      fig <- private$.volcano(
        self$contrast_df,
        self$volcano_spec,
        colour = colour,
        legend = legend,
        scales = scales,
        min_score = min_score
      )
      return(fig)
    },
    #' @description
    #' plotly volcano plots
    #' @param colour column in contrast matrix with colour coding
    #' @return list of ggplots
    #' @param legend default TRUE
    #' @param scales default fixed \code{\link[ggplot2]{facet_wrap}}, scales argument
    #' @param min_score replace p.values or FDR's smaller then min_score with min_score (default 0.0001).
    volcano_plotly = function(
      colour,
      legend = TRUE,
      scales = c("fixed", "free", "free_x", "free_y"),
      min_score = 0.0001
    ) {
      if (missing(colour)) {
        colour <- self$model_name
      }
      scales <- match.arg(scales)
      res <- private$.volcano(
        self$contrast_df,
        self$volcano_spec,
        colour = colour,
        legend = legend,
        scales = scales,
        plotly = TRUE,
        min_score = min_score
      )
      return(res)
    },
    #' @description
    #' ma plot
    #'
    #' MA plot displays the effect size estimate as a function
    #' of the mean protein intensity across groups.
    #' Each dot represents an observed protein.
    #' Red horizontal lines represent the fold-change threshold.
    #'
    #' Sometimes measured effects sizes (differences between samples groups)
    #' are biased by the signal intensity (here protein abundance).
    #' Such systematic effects can be explored using MA-plots.
    #'
    #' @param fc fold change abline
    #' @param colour column in contrast matrix with colour coding
    #' @param legend enable legend default TRUE
    #' @param rank default FALSE, if TRUE then rank of avgAbd is used.
    #' @return ggplot
    ma_plot = function(fc, colour, legend = TRUE, rank = TRUE) {
      if (missing(fc)) {
        fc <- self$fcthresh
      }
      if (missing(colour)) {
        colour <- self$model_name
      }
      private$.ma_fig(self$contrast_df, fc, colour, legend, rank)
    },
    #' @description
    #' ma plotly
    #' @param fc horizontal lines
    #' @param colour column in contrast matrix with colour coding.
    #' @param legend enable legend default TRUE
    #' @param rank default FALSE, if TRUE then rank of avgAbd is used.
    #' @return list of ggplots
    ma_plotly = function(fc, colour, legend = TRUE, rank = FALSE) {
      if (missing(fc)) {
        fc <- self$fcthresh
      }
      if (missing(colour)) {
        colour <- self$model_name
      }
      fig <- private$.ma_fig(self$contrast_df, fc, colour, legend, rank, plotly_mode = TRUE)
      if (!is.null(fig)) {
        fig <- fig |> plotly::ggplotly(tooltip = "subject_id")
      }
      fig
    },
    #' @description
    #' plot a score against the log2 fc e.g. t-statistic
    #' @param scorespec list(score="statistics", fcthres = 2, thresh = 5)
    #' @param colour column with colour coding
    #' @param legend enable legend default TRUE
    #' @return list of ggplots
    score_plot = function(scorespec, colour, legend = TRUE) {
      if (!missing(scorespec)) {
        self$score_spec[[scorespec$score]] <- scorespec
      }
      if (missing(colour)) {
        colour <- self$model_name
      }
      res <- list()
      if (length(self$score_spec) > 0) {
        res <- private$.score_plot(
          self$contrast_df,
          self$score_spec,
          colour = colour,
          legend = legend
        )
      }
      return(res)
    },
    #' @description
    #' plot a score against the log2 fc e.g. t-statistic
    #' @param scorespec list(score="statistics", fcthres = 2, thresh = 5)
    #' @param colour column with colour coding
    #' @param legend enable legend default TRUE
    #' @return list of ggplots
    score_plotly = function(scorespec, colour, legend = TRUE) {
      if (!missing(scorespec)) {
        self$score_spec[[scorespec$score]] <- scorespec
      }
      if (missing(colour)) {
        colour <- self$model_name
      }
      contrast_df <- self$contrast_df |> plotly::highlight_key(~subject_id)
      res <- private$.score_plot(
        contrast_df,
        self$score_spec,
        colour = colour,
        legend = legend
      )

      for (i in seq_along(res)) {
        res[[i]] <- plotly::ggplotly(res[[i]], tooltip = "subject_id")
      }
      return(res)
    },
    #' @description
    #' shows the number of significant proteins per contrasts
    #' @return list which contains ggplots and summary tables
    barplot_threshold = function() {
      bar_results <- list()
      for (i in seq_along(self$volcano_spec)) {
        score_name <- self$volcano_spec[[i]]$score
        score_threshold <- self$volcano_spec[[i]]$thresh
        filt <- dplyr::filter(
          self$contrast_df,
          !is.na(!!sym(score_name)) & !!sym(score_name) < score_threshold
        )
        if (is.numeric(self$fcthresh)) {
          filt <- dplyr::filter(filt, abs(!!sym(self$diff)) > self$fcthresh)
        }
        summary_counts <- group_by(filt, !!sym(self$contrast), !!sym(self$model_name)) |>
          dplyr::summarize(n = n())
        p <- ggplot(summary_counts, aes(x = !!sym(self$contrast), y = n, fill = !!sym(self$model_name))) +
          geom_bar(position = "stack", stat = "identity")
        bar_results[[score_name]] <- list(plot = p, summary = summary_counts)
      }
      return(bar_results)
    }
  ),
  private = list(
    .ma_fig = function(contrast_df, fc, colour, legend, rank, plotly_mode = FALSE) {
      if (is.null(contrast_df[[self$avg.abundance]])) {
        if (!plotly_mode) {
          warning("no group_1 group_2 columns can't generate MA")
        }
        return(NULL)
      }
      abundance_col <- self$avg.abundance
      if (rank) {
        abundance_col <- paste0("rank_", self$avg.abundance)
        contrast_df <- contrast_df |>
          dplyr::group_by(!!sym(self$contrast)) |>
          dplyr::mutate(!!abundance_col := rank(!!sym(self$avg.abundance)))
      } else if (plotly_mode) {
        contrast_df <- contrast_df |> plotly::highlight_key(~subject_id)
      }
      private$.ma_plot(contrast_df, abundance_col, self$diff, self$contrast, fc, colour = colour, legend = legend)
    },
    .volcano = function(
      contrasts,
      scores,
      colour = NULL,
      legend = TRUE,
      scales = "free_y",
      plotly = FALSE,
      min_score = 0.0001
    ) {
      fig <- list()
      for (score in scores) {
        column <- score$score
        contrasts2 <- contrasts |>
          dplyr::filter(!is.na(!!sym(self$diff))) |>
          dplyr::filter(!is.na(!!sym(column))) |>
          dplyr::mutate(!!column := case_when(!!sym(column) < min_score ~ min_score, TRUE ~ !!sym(column)))
        if (plotly) {
          contrasts2 <- contrasts2 |> plotly::highlight_key(~subject_id)
        }
        p <- prolfqua:::.multigroup_volcano(
          contrasts2,
          effect = self$diff,
          significance = column,
          contrast = self$contrast,
          text = "subject_id",
          xintercept = if (is.numeric(self$fcthresh)) {
            c(-self$fcthresh, self$fcthresh)
          } else {
            NULL
          },
          yintercept = score$thresh,
          colour = colour,
          scales = scales
        )
        if (!legend) {
          p <- p + guides(colour = "none")
        }
        if (plotly) {
          p <- plotly::ggplotly(p, tooltip = "subject_id")
        }
        fig[[column]] <- p
      }
      return(fig)
    },
    .histogram = function(score) {
      xlim <- score$xlim
      score <- score$score
      plot <- self$contrast_df |>
        ggplot(aes(x = !!sym(score))) +
        geom_histogram(breaks = seq(from = xlim[1], to = xlim[2], by = xlim[3])) +
        facet_wrap(vars(!!sym(self$contrast)))
      return(plot)
    },
    .ma_plot = function(x, avg.abundance, diff, contrast, fc, colour = NULL, legend = TRUE) {
      p <- ggplot(
        x,
        aes(x = !!sym(avg.abundance), y = !!sym(diff), text = !!sym("subject_id"), colour = !!sym(colour))
      ) +
        geom_point(alpha = 0.5) +
        scale_colour_manual(values = c("black", "green")) +
        facet_wrap(vars(!!sym(contrast)))
      if (FALSE) {
        ylab("log fold change (M)") + xlab("mean log intensities (A)")
      } else {
        NULL
      }
      if (is.numeric(fc)) {
        p <- p + geom_hline(yintercept = c(-fc, fc), linetype = "dashed", colour = "red")
      }

      if (!legend) {
        p <- p + guides(colour = "none")
      }
      return(p)
    },
    .score_plot = function(x, scores, colour, legend = TRUE) {
      fig <- list()
      for (score in scores) {
        xlim <- self$fcthresh
        ylim <- score$thresh
        score <- score$score
        score_values <- if ("data.frame" %in% class(x)) {
          x[[score]]
        } else {
          x$data()[[score]]
        }

        ylims <- c(sign(min(score_values, na.rm = TRUE)) * ylim, sign(max(score_values, na.rm = TRUE)) * ylim)
        p <- ggplot(
          x,
          aes(x = !!sym(self$diff), y = !!sym(score), text = !!sym("subject_id"), colour = !!sym(colour))
        ) +
          scale_colour_manual(values = c("black", "green")) +
          geom_point(alpha = 0.5) +
          facet_wrap(vars(!!sym(self$contrast))) +
          geom_hline(yintercept = c(0), colour = 1) +
          geom_vline(xintercept = c(0), colour = 1) +
          geom_hline(yintercept = ylims, colour = 2, linetype = "dashed")

        if (is.numeric(xlim)) {
          p <- p + geom_vline(xintercept = c(-xlim, xlim), colour = 2, linetype = "dashed")
        }

        if (!legend) {
          p <- p + guides(colour = "none")
        }
        fig[[score]] <- p
      }
      return(fig)
    }
  )
)
